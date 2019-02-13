#!/usr/bin/env python3

"""
Cette version est la version de modification des arbres initialement développées par MM.
Modifiée par G.Louvel pour être importable comme module.

Corrige les arbres d'Ensembl en fonction du seuil minimal de duplication_score et de l'arbre des especes desire
    1: score par defaut (0 -> 1)
    2: coef multiplicateur d'un score reference 6X_species / all_species
    3: duplication_confidence_score calcule sur uniquement 6X_species
"""

from __future__ import print_function

import sys
import itertools

if sys.version_info[0] == 3:
    from LibsDyogen import myFile, myTools, myPhylTree, myProteinTree
else:
    from LibsDyogen.utils import myFile, myTools, myPhylTree, myProteinTree

import logging
logger = logging.getLogger(__name__)


sys.setrecursionlimit(20000)


def alwaysTrue(tree, rnode):
    return True

def alwaysFalse(tree, rnode):
    return False

def setupScoring(phylTree, scoreMethod=1, cutoff=-1):
    """Return a `hasLowScore` function that attributes a return True/False
    depending on whether a duplication node has a good confidence, with
    reference to a given species phylogenetic tree (phyltreefile)."""

    # Limites automatiques de score de duplication
    if scoreMethod in [1, 3]:
        def calc(anc, val):
            return val
    elif scoreMethod == 2:
        def calc(anc, val):
            nesp = len(phylTree.species[anc])
            n2X = len(phylTree.lstEsp2X.intersection(phylTree.species[anc]))
            # La moitie des especes non 2X a vu la duplication (au minimum 1 espece)
            return round(max(1., val*(nesp-n2X)) / nesp, 3) - 2e-3

    minDuplicationScore = {}
    try:
        # Une limite pour tout le monde
        val = float(cutoff)

        # Shortcut
        if val < 0:
            return alwaysFalse

        for anc in phylTree.listAncestr:
            minDuplicationScore[anc] = calc(anc, val)
    except ValueError:
        f = myFile.openFile(cutoff, "r")
        for l in f:
            t = l.split()
            anc = phylTree.officialName[t[0]]
            minDuplicationScore[anc] = calc(anc, float(t[1]))
        f.close()
    logger.debug("minDuplicationScore:\n%s", minDuplicationScore)

    # Les scores dans l'arbre pour les especes modernes valent toujours 1, on
    # doit toujours les accepter
    for esp in phylTree.listSpecies:
        minDuplicationScore[esp] = 0

    @myTools.memoize
    def goodSpecies(anc):
        return phylTree.species[anc].difference(phylTree.lstEsp2X)

    def hasLowScore(tree, rnode):

        logger.debug("# hasLowScore is used.")
        @myTools.memoize
        def getSpeciesSets(node):
            if node in tree.data:
                return set().union(*(getSpeciesSets(x) for (x,_) in tree.data[node]))
            else:
                logger.debug('Node without data (leaf) at %r', tree.info[node]["taxon_name"])
                assert tree.info[node]["taxon_name"] in phylTree.listSpecies
                return set([tree.info[node]["taxon_name"]])

        if rnode not in tree.data:
            return False

        speciessets = [getSpeciesSets(x) for (x,_) in tree.data[rnode]]
        inters = set()
        for (s1,s2) in itertools.combinations(speciessets, 2):
            inters.update(s1.intersection(s2))
        all = set().union(*speciessets)
        anc = tree.info[rnode]["taxon_name"]

        if scoreMethod == 3:
            inters.intersection_update(goodSpecies(anc))
            all.intersection_update(goodSpecies(anc))
        return ((len(inters) == 0) and (minDuplicationScore[anc] == 0)) or (len(inters) < (minDuplicationScore[anc] * len(all)))

    return hasLowScore


def markLowDup(tree, hasLowScore):
    """Mark duplication nodes that have a low score with `'Duplication': 1`.
    Return the count of each type of duplication (dubious, too low, good)."""
    dubious, toolow, good = 0, 0, 0

    # On trie les bonnes duplications des mauvaises
    ################################################
    for (node,inf) in tree.info.items():
        logger.debug('Processing %s: %s', node, inf)
        if inf['Duplication'] != 0:

            if 'dubious_duplication' in inf:
                # On considere que les duplications 'dubious' ne sont pas
                # valables pour une duplication
                assert inf['Duplication'] == 1
                del inf['dubious_duplication']
                dubious += 1

            if hasLowScore(tree, node):
                inf['Duplication'] = 1
                toolow += 1
            else:
                # Il faut la passer a 2 si le score est suffisant
                # Attention: pour les arbres d'Ensembl dont la racine est une
                # duplication, celle-ci vaut 1 (parce qu'elle n'a pas
                # d'outgroup)
                if inf['Duplication'] == 1:
                    inf['Duplication'] = 3
                else:
                    assert inf['Duplication'] in [2,3]
                good += 1
    return dubious, toolow, good


def process(prottrees, phylTree, hasLowScore, defaultFamName="FAM%08d",
            flatten=False, rebuild=False, recurs=True, indicator=False):

    nbEdit = {"dubious": 0, "toolow": 0, "good": 0}
    nbFlattened = 0
    nbRebuilt = 0
    
    nb = 0
    for tree in prottrees:

        assert max(tree.info) < myProteinTree.nextNodeID

        #TODO: only used by the extractGeneFamilies step. Should be ignored here.
        try:
            tree_name = tree.info[tree.root]["tree_name"]
        except KeyError:
            nb += 1
            tree.info[tree.root]["tree_name"] = tree_name = defaultFamName % nb

        dubious, toolow, good = markLowDup(tree, hasLowScore)
        nbEdit["dubious"] += dubious
        nbEdit["toolow"] += toolow
        nbEdit["good"] += good
        
        tree.doBackup()

        if flatten:
            #logger.debug("### Tree %d: First call to flatten ###", tree.backRoot)
            flattened = tree.flattenTree(phylTree, True, indicator=indicator)
            if flattened:
                logger.info("Flattened %5s nodes in tree %6d %r", '?',
                            tree.backRoot, tree_name)
            nbFlattened += int(flattened)

            if tree.backRoot != tree.root:
                logger.warning("Flatten changed root or %r: %d -> %d",
                               tree_name, tree.backRoot, tree.root)

            tree.doBackup()

        if rebuild:
            #assert recurs is False
            
            rebuilt = tree.rebuildTree(phylTree,
                                       hasLowScore if recurs else myProteinTree.alwaysTrue,
                                       indicator=indicator)
            if rebuilt:
                logger.info("Rebuilt   %5s nodes in tree %6d %r", '?',
                            tree.backRoot, tree_name)
            nbRebuilt += int(rebuilt)

            if tree.backRoot != tree.root:
                logger.warning("Rebuild changed root of %r: %d -> %d",
                               tree_name, tree.backRoot, tree.root)
        
        assert "tree_name" in tree.info[tree.root]

        yield tree

    print('%s; %d flattened trees; %d rebuilt trees.' % (
            nbEdit, nbFlattened, nbRebuilt), file=sys.stderr)


if __name__ == '__main__':
    logging.basicConfig()
    logger.setLevel(logging.INFO)

    arguments = myTools.checkArgs(
        [("phylTree.conf",myTools.File), ("ensemblTree",myTools.File)],
        [("flatten",bool,False), ("rebuild",bool,False), ("fam",bool,False),
         ("cutoff",str,"-1"), ("defaultFamName",str,"FAM%08d"),
         ("scoreMethod",int,[1,2,3]), ("newNodeID",int,100000000),
         ("recurs",bool,False), ("indicator",bool,False), ("debug",bool,False)],
        __doc__)
    if arguments['debug']: logger.setLevel(logging.DEBUG)

    myProteinTree.nextNodeID = arguments["newNodeID"]  # For the rebuild step.
    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

    hasLowScore = setupScoring(phylTree,
                               arguments["scoreMethod"],
                               arguments["cutoff"])

    prottrees = myProteinTree.loadTree(arguments["ensemblTree"])

    prottrees = process(prottrees, phylTree, hasLowScore,
                        arguments["defaultFamName"], arguments["flatten"],
                        arguments["rebuild"], arguments["recurs"],
                        arguments["indicator"])

    if arguments["fam"]:
        # Will not work on previous versions of ToolsDyogen.
        from ToolsDyogen.treeTools.ALL.extractGeneFamilies import extractGeneFamilies
        count, dupCount, geneFamilies = extractGeneFamilies(phylTree, prottrees)
    else:
        for tree in prottrees:
            tree.printTree(sys.stdout)

