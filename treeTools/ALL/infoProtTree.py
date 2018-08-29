#! /usr/bin/env python

"""
	Gives information on forest of gene trees.
	Nb families, Nb nodes (speciation, duplication)
"""

from LibsDyogen import myTools, myProteinTree

arguments = myTools.checkArgs([("iniTree", myTools.File)], [], __doc__)

for tree in myProteinTree.loadTree(arguments["iniTree"]):
    next
