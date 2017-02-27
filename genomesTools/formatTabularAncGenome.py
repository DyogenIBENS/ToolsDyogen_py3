#!/usr/bin/env python
# coding: latin-1

import sys
import collections

import utils.myFile
import utils.myPhylTree
import utils.myGenomes
import utils.myTools

__doc__ = """
        transform an ancestral genome in tabular format with descendant species genes (one column by species), with modern position or not

        usage:
                ./formatTabularAncGenome.py PhylTree.conf genome.Boreoeutheria.list.bz2 Boreoeutheria  -in:genesFiles=genes/genesST.%s.list.bz2 +withPos > genome.Boreoeutheria.WithDescendant.list
"""

arguments = utils.myTools.checkArgs(
    [("phylTree.conf", file), ("ancGenome", file), ("target", str)], \
    [("in:genesFiles", str, ""), ("withPos", bool, False)],
    __doc__
)

# loading species tree
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
# Extant species list to load
listSpecies = phylTree.getTargetsSpec(arguments["target"])
newlistSpecies = sorted(listSpecies)

print(utils.myFile.myTSV.printLine(
    ["Anc_chr", "Begin", "End", "Strand", "AncGene", '\t'.join(x for x in newlistSpecies)]), file=sys.stdout)

ancGenome = utils.myGenomes.Genome(arguments["ancGenome"])

genome = {}
for esp in listSpecies:
    # loading extant genome
    if phylTree.isChildOf(esp, arguments["target"]):
        genome[esp] = utils.myGenomes.Genome(arguments["in:genesFiles"] % phylTree.fileName[esp])

desc = {}
for genes in ancGenome:
    strModern = collections.defaultdict(list)
    ancGene1 = genes.names
    #print >> sys.stderr, ancGene1
    desc[ancGene1[0]] = ancGene1[1:]
    #print >> sys.stderr, desc[ancGene1[0]]

    for esp in newlistSpecies:
        strModern[esp]=[]
        for gene in desc[ancGene1[0]]:

            if gene in genome[esp].dicGenes:
                (c, i) = genome[esp].dicGenes[gene]
                # print >>sys.stderr, gene, esp, genome[esp].lstGenes[c][i].chromosome,genome[esp].lstGenes[c][i].beginning,genome[esp].lstGenes[c][i].end
                if arguments["withPos"]:
                    str1 = gene + "/" + str(genome[esp].lstGenes[c][i].chromosome) + ":" + str(
                        genome[esp].lstGenes[c][i].beginning) + "-" + str(genome[esp].lstGenes[c][i].end)
                else:
                    str1 = gene
                strModern[esp].append(str1)

    #print >> sys.stderr, strModern
#     print >> sys.stdout, utils.myFile.myTSV.printLine([genes[0], genes[1], genes[2], genes[3], ancGene1[0], "\t".join(
 #       str(x[1][1:])[1:-1] for x in sorted(strModern.iteritems(), reverse=False))])
    print(utils.myFile.myTSV.printLine([genes[0], genes[1], genes[2], genes[3], ancGene1[0], "\t".join(
        str(x[1])[1:-1] for x in sorted(iter(strModern.items()), reverse=False))]), file=sys.stdout)
