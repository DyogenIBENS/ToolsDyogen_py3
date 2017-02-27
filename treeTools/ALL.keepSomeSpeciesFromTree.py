#! /usr/bin/env python

__doc__ = """
	Load the forest of gene trees and reformat it in order to keep only the requested species
		Usage:
			./ALL.keepSomeSpeciesFromTree.py PhylTree.conf GeneTrees.bz2 Species.txt > Subtrees
"""

import os
import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs([("phylTree.conf", file), ("proteinTree", file), ("SpeciesList", file)], [], __doc__)

count = collections.defaultdict(list)
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
allTrees = {}
allRoots = []

for tree in utils.myProteinTree.loadTree(arguments["proteinTree"]):

    allTrees[tree.root] = tree
    allRoots.append(tree.root)


    def countGenes(node):

        if node in tree.data:
            for (g, _) in tree.data[node]:
                countGenes(g)
        else:
            count[(tree.info[node]['gene_name'], tree.info[node]['taxon_name'])].append((node, tree.root))


    countGenes(tree.root)

SpeciesList1 = []
for line in open(arguments["SpeciesList"]):
    ch = line.split('\n')
    SpeciesList1.append(ch[0])

SpeciesList = []
for esp in phylTree.listSpecies:
    # print 'esp', esp
    if esp not in SpeciesList1:
        SpeciesList.append(esp)

print('keep', SpeciesList1, file=sys.stderr)
print('remove', SpeciesList, file=sys.stderr)

for ((name, esp), l) in count.items():
    # print >> sys.stderr, name, esp, l
    if esp not in SpeciesList:
        continue
    delete = True
    s = set(root for (node, root) in l)
    # print >> sys.stderr, name, esp, len(l), len(s), s
    if len(s) == 1:
        root = s.pop()
        lca = None
        tree = allTrees[root]


        def findLastCommonAncestor(node):
            global lca
            global delete
            if node in tree.data:
                nOK = nTot = 0
                for (g, _) in tree.data[node]:
                    res = findLastCommonAncestor(g)
                    nOK += res[0]
                    nTot += res[1]
                # if (nOK == nTot) and (nOK == len(l)):
                #	lca = node
                #	delete = False
                return (nOK, nTot)

            else:
                if tree.info[node]['gene_name'] == name:
                    return (1, 1)
                else:
                    return (0, 1)


        findLastCommonAncestor(root)
    # print >> sys.stderr, delete
    # Rebuild the tree
    if delete:
        for (node, root) in l:
            tree = allTrees[root]


            def recdelete(x):
                if x in tree.data:
                    ll = [(g, l) for (g, l) in tree.data[x] if recdelete(g)]
                    if len(ll) == 0:
                        del tree.data[x]
                    else:
                        tree.data[x] = ll
                    return x in tree.data
                # data[x] = [(g,l) for (g,l) in data[x] if g != node]
                # for (g,_) in data[x]:
                #	recdelete(g)
                else:
                    return x != node


            recdelete(root)
    else:
        del tree.data[lca]
        tree.info[lca] = {
            ['taxon_name']: esp,
            ['Duplication']: 0,
            ['gene_name']: name,
            ['transcript_name']: "/".join(info[node]['transcript_name'] for (node, _) in l),
            ['protein_name']: "/".join(info[node]['protein_name'] for (node, _) in l)
        }

for root in allRoots:
    tree = allTrees[root]

    if tree.info[root]['taxon_name'] not in SpeciesList:
        # if gene == 1:
        tree.flattenTree(phylTree, True)

allRoots2 = []
for i in allRoots:
    # print >> sys.stderr, 'i', i
    if allTrees[i].data:
        allRoots2.append(i)

for root in allRoots2:
    tree = allTrees[root]
    tree.flattenTree(phylTree, True)
    tree.printTree(sys.stdout)
