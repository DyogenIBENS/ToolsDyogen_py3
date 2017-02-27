#! /usr/bin/env python

__doc__ = """
	Decoupe un fichier d'arbres en fichiers separes

	usage:
			./splitTrees.py GeneTreeForest.phylTree.bz2 Fam.%s
"""

import utils.myFile
import utils.myTools
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("proteinTree",file), ("output",str)], [], __doc__ )

for (i,tree) in enumerate(utils.myProteinTree.loadTree(arguments["proteinTree"])):
	print(i)
	f = utils.myFile.openFile(arguments["output"] % (i+1), "w")
	tree.printTree(f)
	f.close()

