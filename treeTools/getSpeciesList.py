#! /usr/bin/env python

__doc__ = """
	Give all the descendant extant species for a given ancestor in a species tree
	
	Usage:	getSpeciesList.py PhylTree.conf Boreoeutheria
"""

import utils.myFile
import utils.myTools
import utils.myPhylTree

# Arguments
arguments = utils.myTools.checkArgs( \
    [("phylTree.conf", file), ("anc", str)], [], \
    __doc__ \
    )


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for (x, _) in phylTree.items[arguments["anc"]]:
    for y in phylTree.species[x]:
        print(y)
