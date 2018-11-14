#! /usr/bin/env python3

"""
	Extract Newick or NHX trees from Phyltree protein trees

	usage:
		./ALL.extractNewickTrees.py GeneTrees.bz2 +withDist +withNHXTags +withAncSpeciesNames +withAncGenesNames
"""

# Librairies
import sys

import LibsDyogen.myTools       as myTools
import LibsDyogen.myProteinTree as myProteinTree

# Arguments
arguments = myTools.checkArgs( [("proteinTree",myTools.File)],
                               [("withDist",bool,False),
                                ("withNHXTags",bool,False),
                                ("withAncSpeciesNames",bool,False),
                                ("withAncGenesNames",bool,False)], __doc__ )

print("Mise en forme des arbres ...", end=' ', file=sys.stderr)
nb = 0
for tree in myProteinTree.loadTree(arguments["proteinTree"]):
	tree.printNewick(sys.stdout, withDist=arguments["withDist"], withTags=arguments["withNHXTags"], withAncSpeciesNames=arguments["withAncSpeciesNames"], withAncGenesNames=arguments["withAncGenesNames"])
	nb += 1
print("%d arbres OK" % nb, file=sys.stderr)

