#!/usr/bin/env python3

"""
	Convertit un genome (suite de diagonales) en genome (suite de genes)
"""

import sys

import itertools
from LibsDyogen import myTools, myGenomes

arguments = myTools.checkArgs( [("contigsFile",file), ("ancGenesFile",file)], [], __doc__)

ancGenes = myGenomes.Genome(arguments["ancGenesFile"])

genome = myGenomes.Genome(arguments["contigsFile"], ancGenes=ancGenes)

genome.printEnsembl(sys.stdout)

