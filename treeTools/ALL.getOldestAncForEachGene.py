#! /usr/bin/env python

__doc__ = """
	Print the oldest ancestor of each request species gene.
	
	usage:
		./ALL.getOldestAncForEachGene.py Foret_arbre_ensembl82_0.30.bz2 "Homo sapiens"
"""

import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs([("proteinTree", file), ("taxon_name", str)], [], __doc__)



# Information on Gene
def printGeneNode(node):
    txt = [tree.info[tree.root]["taxon_name"]]
    #txt.append("GENE")
    #txt.append(tree.info[node].pop("taxon_name", None))
    txt.append(tree.info[node].pop("gene_name", None))
    print(utils.myFile.myTSV.printLine(txt))


# Recursive loop on the gene family
def do(node):
    if node in tree.data:
        for (g, d) in tree.data[node]:
		do(g)
    elif tree.info[node]["taxon_name"] == arguments["taxon_name"]:
        printGeneNode(node)


# searching for the good gene tree
for tree in utils.myProteinTree.loadTree(arguments["proteinTree"]):
    
    do(tree.root)
