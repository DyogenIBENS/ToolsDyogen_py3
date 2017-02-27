#! /usr/bin/env python

__doc__ = """
	Cut the gene trees under a given ancestral species
	
	usage:
		ALL.cutTree.py PhylTree.conf GeneTrees.bz2 Amniota > Amniota.genetrees
"""

import sys

import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs([("phylTree.conf", file), ("iniTree", file), ("rootSpecies", str)], [], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


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
for tree in utils.myProteinTree.loadTree(arguments["iniTree"]):
    l = search(tree.root)
    nb += len(l)
    if len(l) == 1:
        tree.info[l[0]]["tree_name"] = tree.info[tree.root]["tree_name"]
        utils.myProteinTree.printTree(sys.stdout, tree.data, tree.info, l[0])
    else:
        for (i, r) in enumerate(l):
            tree.info[r]["tree_name"] = tree.info[tree.root]["tree_name"] + utils.myProteinTree.getDupSuffix(i + 1,
                                                                                                             True)
            utils.myProteinTree.printTree(sys.stdout, tree.data, tree.info, r)

print(nb, "extracted trees", file=sys.stderr)
