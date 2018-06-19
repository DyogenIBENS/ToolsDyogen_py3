#!/usr/bin/env python3

from __future__ import print_function

"""
	Extract Trees that contain Specific ancestral nodes and/or do not contain specific ancestral nodes

	usage:
		./ALL.extractSpecificTrees.py GeneTrees.bz2 -inDesc='Sauria,Percomorpha' -notinDesc='Mammalia'
"""

import sys
import collections

import LibsDyogen.myFile        as myFile
import LibsDyogen.myTools       as myTools
import LibsDyogen.myPhylTree    as myPhylTree
import LibsDyogen.myProteinTree as myProteinTree
import LibsDyogen.myGenomes     as myGenomes

arguments = myTools.checkArgs([("phylTree.conf", file),
                                     ("proteinTree", file)],
                                    [("inDesc", str, "Mammalia"),
                                     ("notinDesc", str, "HomoPan"),
                                     ("in:ancGenesFiles", str, "")],
                                    __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


lstIn = arguments["inDesc"].split(",")

print(lstIn, file=sys.stderr)

lstOut = arguments["notinDesc"].split(",")

print(lstOut, file=sys.stderr)

lst = lstIn+lstOut

print(lst, file=sys.stderr)

lastCommonAnc = phylTree.lastCommonAncestor([x for x in lst])

print("LastCommonAnc = ", lastCommonAnc, file=sys.stderr)

# Chargement des fichiers
lastAncGenes = myGenomes.Genome(arguments["in:ancGenesFiles"] % lastCommonAnc)

in_genome  = {}
out_genome = {}
for inAnc in lstIn:
    in_genome[inAnc]=[]
    for gene in myGenomes.Genome(arguments["in:ancGenesFiles"] % inAnc):
        in_genome[inAnc].append(gene.names[0])
    in_genome[inAnc]=set(in_genome[inAnc])

#print in_genome
intersectIn=set.intersection(*(list(in_genome.values())))

#print intersectIn
print(len(intersectIn), file=sys.stderr)


for notinAnc in  lstOut:
    print(notinAnc, file=sys.stderr)
    out_genome[notinAnc]=[]
    for gene in myGenomes.Genome(arguments["in:ancGenesFiles"] % notinAnc):
        out_genome[notinAnc].append(gene.names[0])
    out_genome[notinAnc]=set(out_genome[notinAnc])

intersectOut=set.intersection(*(list(out_genome.values())))
print(len(intersectOut), file=sys.stderr)


for gene in lastAncGenes:
    if (gene.names[0] in intersectIn) and (gene.names[0] not in intersectOut):
        print(gene.names[0], " ".join((x for x in gene.names[1:])))

"""
for ancGene in lastAncGenes:

    if (ancGene.names[0] in in_genome):
        print ancGene.names[0]
"""

"""
# parsing trees
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


allTrees = {}
allRoots = []
# searching for the good gene tree

for tree in myProteinTree.loadTree(arguments["proteinTree"]):
    allTrees[tree.root] = tree
    allRoots.append(tree.root)
    count = collections.defaultdict(list)


    def countGenes(node):

        if node in tree.data:
            #print "node: ", node
            for (g, _) in tree.data[node]:
            #    print "g: ", g
                countGenes(g)
            #print tree.info[node]['taxon_name']
            count[(tree.info[node]['taxon_name'])].append((node, tree.root))



    countGenes(tree.root)
    if ( not( count['Mammalia'] ) and count['Sauria'] and count['Neopterygii']):
        print "found: ", tree.data[tree.root]
        tree.printTree(sys.stdout)

"""
