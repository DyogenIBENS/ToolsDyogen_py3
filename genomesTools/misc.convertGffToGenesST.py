#!/usr/bin/env python
# coding: latin-1

__doc__ = """
    Transform a reconstructed ancestral list of blocks (diags) in a formated ancestral genome (list of genes)

    usage:
        ./misc.convertGffToGenesST.py.py gffFile
"""

import sys

from gff3 import Gff3
import utils.myGenomes
import utils.myTools



arguments = utils.myTools.checkArgs([("gffFile", file)], [], __doc__)

gff = Gff3(arguments["gffFile"])
genes = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'mRNA']

for gene in genes:
    #print >> sys.stdout,  gene['seqid'], gene['start'], gene['end'], gene['strand'], gene['attributes']['ID']
    if gene['strand'] == "+":
        gene['strand'] = '1'
    else:
        gene['strand'] = '-1'
    print(utils.myFile.myTSV.printLine([gene['seqid'], gene['start'], gene['end'], gene['strand'], gene['attributes']['ID']]), file=sys.stdout)