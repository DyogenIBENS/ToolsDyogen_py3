#! /usr/bin/env python

__doc__ = """
	Convertit un genome (suite de diagonales) en genome (suite de genes)
"""

import sys

import itertools
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [("contigsFile",file), ("ancGenesFile",file)], [], __doc__)

ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"])

genome = utils.myGenomes.Genome(arguments["contigsFile"], ancGenes=ancGenes)

genome.printEnsembl(sys.stdout)

