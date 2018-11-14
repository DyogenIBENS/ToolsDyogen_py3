#! /usr/bin/env python3

"""
	Print the history of a gene giving its gene tree.
	
	usage:
		./ALL.getGeneHistory.py GeneTrees.bz2 ENSOARG00000011785
"""

import sys

from LibsDyogen import myFile, myTools, myPhylTree, myProteinTree


arguments = myTools.checkArgs([("proteinTree", myTools.File),
                               ("gene_name", str)],
                              [], __doc__)


# Information on ancestral node
def printAncNode(node):
    txt = [node]
    d = tree.info[node].pop('Duplication', None)
    if tree.info[node].pop("dubious_duplication", None):
        txt.append("DUBIOUS_DUPLICATION")
    elif (d == 1) and ("duplication_confidence_score" in tree.info[node]):
        txt.append("ROOT_DUPLICATION")
    elif d == 2:
        txt.append("DUPLICATION")
    elif d == 3:
        txt.append("EDITED_DUPLICATION")
    else:
        txt.append("SPECIATION")
    txt.append(tree.info[node].pop("taxon_name", None))
    txt.append(tree.info[node].pop("family_name", None))
    txt.append(tree.info[node].pop("Bootstrap", None))
    txt.append(tree.info[node].pop("duplication_confidence_score", None))
    print(myFile.myTSV.printLine(txt))


# Information on Gene
def printGeneNode(node):
    txt = [node]
    txt.append("GENE")
    txt.append(tree.info[node].pop("taxon_name", None))
    txt.append(tree.info[node].pop("gene_name", None))
    print(myFile.myTSV.printLine(txt))


# Recursive loop on the gene family
def do(node):
    if node in tree.data:
        for (g, d) in tree.data[node]:
            if do(g):
                printAncNode(node)
                return True
    elif tree.info[node]["gene_name"] == arguments["gene_name"]:
        printGeneNode(node)
        return True
    return False


# searching for the good gene tree
for tree in myProteinTree.loadTree(arguments["proteinTree"]):
    if do(tree.root):
        break
