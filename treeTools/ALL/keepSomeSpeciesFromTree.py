#!/usr/bin/env python3


"""
Load the forest of gene trees and reformat it in order to keep only the requested species

Usage:
    ./ALL.keepSomeSpeciesFromTree.py PhylTree.conf GeneTrees.bz2 Species.txt > Subtrees

Options:
    +invert: remove the species from the given list.
"""


import os
import sys
import collections

from LibsDyogen import myTools, myPhylTree, myProteinTree


arguments = myTools.checkArgs([("phylTree.conf", myTools.File),
                               ("proteinTree",   myTools.File),
                               ("SpeciesList",   myTools.File)],
                              [("invert", bool, False)],
                              __doc__)

count = collections.defaultdict(list)
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
allTrees = {}
allRoots = []


#def countGenes(tree, node):
#    if node in tree.data:
#        for (g, _) in tree.data[node]:
#            countGenes(tree, g)
#    else:
#        count[(tree.info[node]['gene_name'],
#               tree.info[node]['taxon_name'])].append((node, tree.root))

def countGenes(tree):
    for leaf in set(tree.info).difference(tree.data):
        try:
            count[(tree.info[leaf]['gene_name'],
                   tree.info[leaf]['taxon_name'])].append((leaf, tree.root))
        except KeyError as err:
            err.args += (tree.root, leaf)
            raise


for tree in myProteinTree.loadTree(arguments["proteinTree"]):
    allTrees[tree.root] = tree
    allRoots.append(tree.root)
    countGenes(tree)


with open(arguments["SpeciesList"]) as species_f:
    keptSpecies = [line.rstrip() for line in species_f]


if arguments["invert"]:
    removedSpecies = keptSpecies
    keptSpecies = phylTree.listSpecies.difference(removedSpecies)
else:
    removedSpecies = phylTree.listSpecies.difference(keptSpecies)


print('INFO:Keep', ' '.join(keptSpecies), file=sys.stderr)
print('INFO:Remove', ' '.join(removedSpecies), file=sys.stderr)


# Iterate over all bad leaves (not in order)
for ((name, esp), l) in count.items():
    
    if esp not in removedSpecies:
        continue  #So why putting it in `count` in the first place?

    delete = True
    s = set(root for (node, root) in l)
    # print >> sys.stderr, name, esp, len(l), len(s), s
    assert len(s) == 1, 'Expecting 1 root for a given (gene_name, taxon_name)'
    assert len(l) == 1, 'Expecting 1 (root,leaf) for a given (gene_name, taxon_name)'

    root = s.pop()
    lca = None
    tree = allTrees[root]

    # Unused.
    def findLastCommonAncestor(node):
        # It its present state, does nothing but count.
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
    assert delete, "'delete' Should not have changed"
    if delete:
        for (badleaf, root) in l:
            tree = allTrees[root]

            def recdelete(x):
                """Find the bad child (badleaf), and signal it to its parent."""
                if x in tree.data:

                    # List only the good children (not supporting the badleaf):
                    ll = [(g, l) for (g, l) in tree.data[x] if recdelete(g)]
                    
                    if len(ll) == 0:
                        del tree.data[x]
                    else:
                        tree.data[x] = ll

                    return x in tree.data
                # data[x] = [(g,l) for (g,l) in data[x] if g != badleaf]
                # for (g,_) in data[x]:
                #	recdelete(g)
                else:
                    return x != badleaf

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

    if tree.info[root]['taxon_name'] not in removedSpecies:
        # if gene == 1:
        tree.flattenTree(phylTree, True)

allRoots2 = []
for i in allRoots:
    # print >> sys.stderr, 'i', i
    if allTrees[i].data:
        allRoots2.append(i)

for root in allRoots2:
    tree = allTrees[root]
    # PB: this recdelete+flatten doesn't get the new branch length correctly:
    # assigns deleted branch length to children items.
    tree.flattenTree(phylTree, True)
    tree.printTree(sys.stdout)
