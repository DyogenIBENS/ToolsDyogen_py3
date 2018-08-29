#! /usr/bin/env python

"""
	Decoupe un fichier d'arbres en fichiers separes

	usage:
			./splitTrees.py GeneTreeForest.phylTree.bz2 Fam.%s
"""

from LibsDyogen import myFile, myTools, myProteinTree

arguments = myTools.checkArgs([("proteinTree", myTools.File), ("output",str)], 
                              [], __doc__ )

for (i,tree) in enumerate(myProteinTree.loadTree(arguments["proteinTree"])):
	print(i)
	f = myFile.openFile(arguments["output"] % (i+1), "w")
	tree.printTree(f)
	f.close()

