#! /usr/bin/env python3

"""
    Give all the descendant extant species for a given ancestor in a species tree
    
    Usage: getSpeciesList.py PhylTree.conf [Boreoeutheria]
"""

import LibsDyogen.myTools    as myTools
import LibsDyogen.myPhylTree as myPhylTree

# Arguments
arguments = myTools.checkArgs( \
    [("phylTree.conf", myTools.File)], [("anc", str, 'root')], \
    __doc__ \
    )


# L'arbre phylogenetique
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

anc = phylTree.root if arguments["anc"] == 'root' else arguments["anc"]

for sp in sorted(phylTree.species[anc]):
    print(sp)
