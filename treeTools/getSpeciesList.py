#! /usr/bin/env python3

"""
    Give all the descendant extant species for a given ancestor in a species tree
    
    Usage: getSpeciesList.py PhylTree.conf Boreoeutheria
"""

import LibsDyogen.myTools    as myTools
import LibsDyogen.myPhylTree as myPhylTree

# Arguments
arguments = myTools.checkArgs( \
    [("phylTree.conf", myTools.File), ("anc", str)], [], \
    __doc__ \
    )


# L'arbre phylogenetique
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for (x, _) in phylTree.items[arguments["anc"]]:
    for y in phylTree.species[x]:
        print(y)
