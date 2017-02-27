#! /usr/bin/env python

__doc__ = """
        Count of genes per chromosome in a genome file
	
	Usage: 	ALL.statsKaryotype.py genesST.Homo.sapiens.list.bz2 -minChrSize=20
		ALL.statsKaryotype.py genome.Boreoeutheria.list.bz2 
"""

import sys

import utils.myGenomes
import utils.myTools

arguments = utils.myTools.checkArgs([("genesFiles", str)], [("minChrSize", int, 1)], __doc__)

genome = utils.myGenomes.Genome(arguments["genesFiles"])

# print >> sys.stderr, genome
# print >> sys.stdout, "Chr","Length"
for (chrom, l) in genome.lstGenes.items():
    if len(l) >= arguments["minChrSize"]:
        print(chrom, len(l), file=sys.stdout)
