#! /usr/bin/env python3

"""
Extract ancestral gene content from a forest of gene trees.
One file per ancestor, and one file per extant species.

usage:
    ./treeTools/ALL/extractGeneFamilies.py PhylTree.conf GeneTrees.bz2 -out:ancGenesFiles=ancGenes/all/ancGenes.%s.list.bz2
    ./treeTools/ALL/extractGeneFamilies.py PhylTree.conf GeneTrees.bz2 -out:ancGenesFiles=ancGenes/all/ancGenes.%s.list.bz2 +bz2 > ForestWithFamilyNames.bz2
"""

from __future__ import print_function


import collections
import sys

from LibsDyogen import myFile, myPhylTree, myProteinTree, myTools

sys.setrecursionlimit(10000)


### TODO: Put in LibsDyogen
### TODO: Make a subfunction that applies only on one protein tree at once.
def extractGeneFamilies(phylTree, proteinTrees, reuseNames=False, outFile=sys.stdout):

    get_baseName = (lambda baseName: baseName) if reuseNames else \
                   (lambda baseName: baseName.split(".")[0])

    dupCount = collections.defaultdict(int)
    def futureName(name, dup):
        if dup >= 2:
            dupCount[name] += 1
            return name + myProteinTree.getDupSuffix(dupCount[name], False)
        else:
            return name

    def getRoots(node, previousAnc, lastWrittenAnc):
        """Find true root in each family"""
        ### Unused function
        newAnc = tree.info[node]['taxon_name']
        _, newLastWritten, isroot = myProteinTree.getIntermediateAnc(phylTree,
                                        previousAnc, lastWrittenAnc, newAnc,
                                        tree.info[node]['Duplication'] >= 2)

        if isroot:
            return [node]

        # descendant genes
        subRoots = []
        for (g, _) in tree.data.get(node, []):
            subRoots.extend(getRoots(g, newAnc, newLastWritten))
        return subRoots

    count = collections.defaultdict(int)
    geneFamilies = collections.defaultdict(list)

    def doGeneFamilies(node, baseName, previousAnc, lastWrittenAnc):
        """Backup all the genes families"""
        newAnc = tree.info[node]['taxon_name']
        toWrite, newLastWritten, isroot = myProteinTree.getIntermediateAnc(phylTree,
                                                   previousAnc, lastWrittenAnc,
                                                   newAnc, tree.info[node]['Duplication'] >= 2)

        if isroot and (previousAnc is not None):
            baseName = get_baseName(baseName)
            count[baseName] += 1
            currName = baseName + myProteinTree.getDupSuffix(count[baseName], True)
        else:
            currName = baseName
        tree.info[node]['family_name'] = currName

        # descendant genes
        # print >> sys.stderr, tree.info[node]
        if node in tree.data:
            allGenes = []
            for (g, _) in tree.data[node]:
                try:
                    allGenes.extend(
                        doGeneFamilies(g,
                                futureName(currName,
                                           tree.info[node]['Duplication']),
                                newAnc,
                                newLastWritten))
                except BaseException as err:
                    err.args += ("Child id, node id: (%d, %d)" % (g,node),)
                    raise
        else:
            allGenes = [tree.info[node]["gene_name"]]

        for a in toWrite:
            geneFamilies[a].append([currName] + allGenes)

        return allGenes

    for tree in proteinTrees:
        # print >> sys.stderr, tree.root
        try:
            doGeneFamilies(tree.root, tree.info[tree.root]["tree_name"], None, None)
        except BaseException as err:
            err.args += ("Root id: '%d'" % tree.root,)
            raise
        tree.printTree(outFile)

    return count, dupCount, geneFamilies


if __name__ == '__main__':

    # Arguments
    arguments = myTools.checkArgs([("phylTree.conf", myTools.File),
                                   ("proteinTree", myTools.File)],
                                  [("out:ancGenesFiles", str, ""),
                                   ("reuseNames", bool, False)],
                                  __doc__)

    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
    proteinTrees = myProteinTree.loadTree(arguments["proteinTree"])

    count, dupCount, geneFamilies = extractGeneFamilies(phylTree,
                                                        proteinTrees,
                                                        arguments["reuseNames"])

    outTemplate = arguments["out:ancGenesFiles"]
    for (anc, lst) in geneFamilies.items():
        print("Ecriture des familles de %s ..." % anc, end=' ', file=sys.stderr)
        f = myFile.openFile(outTemplate % phylTree.fileName[anc], "w")
        for gg in lst:
            print(" ".join(gg), file=f)
        f.close()
        print(len(lst), "OK", file=sys.stderr)
