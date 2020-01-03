#! /usr/bin/env python

"""
	Cut the gene trees under a given ancestral species
	
	usage:
		ALL.cutTree.py PhylTree.conf GeneTrees.bz2 Amniota > Amniota.genetrees
"""

import sys

from LibsDyogen import myTools, myPhylTree, myProteinTree


def main():
    arguments = myTools.checkArgs([("phylTree.conf", myTools.File),
                                   ("iniTree", myTools.File),
                                   ("rootSpecies", str)], [], __doc__)

    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


    # Returns a list of nodes under the new root species
    #########################################################
    def search(node):
        if phylTree.isChildOf(tree.info[node]['taxon_name'], arguments["rootSpecies"]):
            return [node]
        elif node in tree.data:
            r = []
            for (g, _) in tree.data[node]:
                r.extend(search(g))
            return r
        else:
            return []


    nb = 0
    for tree in myProteinTree.loadTree(arguments["iniTree"]):
        l = search(tree.root)
        nb += len(l)
        if len(l) == 1:
            tree.info[l[0]]["tree_name"] = tree.info[tree.root]["tree_name"]
            myProteinTree.printTree(sys.stdout, tree.data, tree.info, l[0])
        else:
            for (i, r) in enumerate(l):
                tree.info[r]["tree_name"] = tree.info[tree.root]["tree_name"] + myProteinTree.getDupSuffix(i + 1,
                                                                                                                 True)
                myProteinTree.printTree(sys.stdout, tree.data, tree.info, r)

    print(nb, "extracted trees", file=sys.stderr)


if __name__ == '__main__':
    main()
