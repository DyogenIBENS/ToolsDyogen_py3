#!/usr/bin/env python3

"""
        Count of genes per chromosome in a genome file

    Usage:     ALL.statsKaryotype.py genesST.Homo.sapiens.list.bz2 -minChrSize=20
        ALL.statsKaryotype.py genome.Boreoeutheria.list.bz2 
"""

import sys

from LibsDyogen import myGenomes, myTools

arguments = myTools.checkArgs([("genesFiles", str)], [("minChrSize", int, 1)], __doc__)

genome = myGenomes.Genome(arguments["genesFiles"])

# print >> sys.stderr, genome
# print >> sys.stdout, "Chr","Length"
for (chrom, l) in genome.lstGenes.items():
    if len(l) >= arguments["minChrSize"]:
        print(chrom, len(l), file=sys.stdout)
