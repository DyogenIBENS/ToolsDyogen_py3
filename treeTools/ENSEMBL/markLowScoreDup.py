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

sys.setrecursionlimit(20000)


def alwaysTrue(tree, rnode):
    return True


def setup_scoring(phylTree, scoreMethod=1, cutoff=-1):
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
        for anc in phylTree.listAncestr:
            minDuplicationScore[anc] = calc(anc, val)
    except ValueError:
        f = myFile.openFile(cutoff, "r")
        for l in f:
            t = l.split()
            anc = phylTree.officialName[t[0]]
            minDuplicationScore[anc] = calc(anc, float(t[1]))
        f.close()
    print("minDuplicationScore:", minDuplicationScore, file=sys.stderr)

    # Les scores dans l'arbre pour les especes modernes valent toujours 1, on
    # doit toujours les accepter
    for esp in phylTree.listSpecies:
        minDuplicationScore[esp] = 0

    @myTools.memoize
    def goodSpecies(anc):
        return phylTree.species[anc].difference(phylTree.lstEsp2X)

    def hasLowScore(tree, rnode):

        @myTools.memoize
        def getSpeciesSets(node):
            if node in tree.data:
                return set().union(*(getSpeciesSets(x) for (x,_) in tree.data[node]))
            else:
                print(tree.info[node]["taxon_name"], file=sys.stderr)
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


def markLowDup(tree, hasLowScore, newNodeID=100000000):
    """Mark duplication nodes that have a low score with `'Duplication': 1`.
    Return the count of each type of duplication (dubious, too low, good)."""
    dubious, toolow, good = 0, 0, 0
    assert max(tree.info) < newNodeID

    # On trie les bonnes duplications des mauvaises
    ################################################
    for (node,inf) in tree.info.items():
        print(node,inf, file=sys.stderr)
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


if __name__ == '__main__':
    arguments = myTools.checkArgs( \
        [("phylTree.conf",myTools.File), ("ensemblTree",myTools.File)], \
        [("flatten",bool,False), ("rebuild",bool,False),
         ("cutoff",str,"-1"), ("defaultFamName",str,"FAM%08d"),
         ("scoreMethod",int,[1,2,3]), ("newNodeID",int,100000000),
         ("recurs",bool,False)], \
        __doc__ \
    )

    myProteinTree.nextNodeID = arguments["newNodeID"]
    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

    hasLowScore = setup_scoring(phylTree,
                                arguments["scoreMethod"],
                                arguments["cutoff"])

    nbEdit = {"dubious": 0, "toolow": 0, "good": 0}
    nbFlattened = 0
    nbRebuilt = 0
    

    for (nb,tree) in enumerate(myProteinTree.loadTree(arguments["ensemblTree"])):

        dubious, toolow, good = markLowDup(tree, hasLowScore, arguments["newNodeID"])
        nbEdit["dubious"] += dubious
        nbEdit["toolow"] += toolow
        nbEdit["good"] += good
        if arguments["flatten"]:
            flattened = tree.flattenTree(phylTree, True)
            nbFlattened += int(flattened)
        if arguments["rebuild"]:
            rebuilt = tree.rebuildTree(phylTree, hasLowScore if arguments["recurs"] else alwaysTrue)
            nbRebuilt += int(rebuilt)

        if "tree_name" not in tree.info[tree.root]:
            tree.info[tree.root]["tree_name"] = arguments["defaultFamName"] % nb

        tree.printTree(sys.stdout)

    print(nbEdit, 'flattened:', nbFlattened, 'rebuilt:', nbRebuilt, file=sys.stderr)
