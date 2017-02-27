#! /usr/bin/env python

__doc__ = """
	Gives information on forest of gene trees.
	Nb families, Nb nodes (speciation, duplication)
"""

import utils.myTools
import utils.myProteinTree

arguments = utils.myTools.checkArgs([("iniTree", file)], [], __doc__)

for tree in utils.myProteinTree.loadTree(arguments["iniTree"]):
    next
