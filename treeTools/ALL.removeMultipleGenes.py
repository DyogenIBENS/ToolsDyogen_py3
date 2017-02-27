#! /usr/bin/env python

__doc__ = """
	Analyse l'arbre et supprime les genes qui sont presents en plusieurs copies
"""

import os
import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

sys.setrecursionlimit(10000)
arguments = utils.myTools.checkArgs([("phylTree.conf", file), ("proteinTree", file)], [], __doc__)

count = collections.defaultdict(list)

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

for ((name, esp), l) in count.items():
    if len(l) == 1:
        continue
    delete = True
    s = set(root for (node, root) in l)
    print(name, esp, len(l), len(s), end=' ', file=sys.stderr)
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
                if (nOK == nTot) and (nOK == len(l)):
                    lca = node
                    delete = False
                return (nOK, nTot)

            else:
                if tree.info[node]['gene_name'] == name:
                    return (1, 1)
                else:
                    return (0, 1)


        findLastCommonAncestor(root)

    print(delete, file=sys.stderr)
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
            'taxon_name': esp,
            'Duplication': 0,
            'gene_name': name,
            'transcript_name': "/".join(info[node]['transcript_name'] for (node, _) in l),
            'protein_name': "/".join(info[node]['protein_name'] for (node, _) in l)
        }

for root in allRoots:
    allTrees[root].printTree(sys.stdout)
