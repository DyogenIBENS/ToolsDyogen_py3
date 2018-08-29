#! /usr/bin/env python

"""
	Print the oldest ancestor of each request species gene.
	
	usage:
		./ALL.getOldestAncForEachGene.py Foret_arbre_ensembl82_0.30.bz2 "Homo sapiens"
"""

import sys

from LibsDyogen import myFile, myTools, myPhylTree, myProteinTree

arguments = myTools.checkArgs([("proteinTree", myTools.File),
                               ("taxon_name", str)], [], __doc__)



# Information on Gene
def printGeneNode(node):
    txt = [tree.info[tree.root]["taxon_name"]]
    #txt.append("GENE")
    #txt.append(tree.info[node].pop("taxon_name", None))
    txt.append(tree.info[node].pop("gene_name", None))
    print(myFile.myTSV.printLine(txt))


# Recursive loop on the gene family
def do(node):
    if node in tree.data:
        for (g, d) in tree.data[node]:
		do(g)
    elif tree.info[node]["taxon_name"] == arguments["taxon_name"]:
        printGeneNode(node)


# searching for the good gene tree
for tree in myProteinTree.loadTree(arguments["proteinTree"]):
    
    do(tree.root)
