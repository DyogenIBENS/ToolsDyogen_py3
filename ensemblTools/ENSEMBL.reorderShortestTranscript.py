#!/usr/bin/env python3


"""
    Creates a list of ordered genes according to the shorstest transcript for genes with alternatve splicing transcript
    Usage:

        ENSEMBL.reorderShortestTranscript.py orig/genes.Homo.sapiens.list.bz2 coords/coordsHomo.sapiens.list.bz2 -authorizedBiotypes=protein_coding | bzip2 > genesST.Homo.sapiens.list.bz2
"""

import collections
import sys

from LibsDyogen import myFile, myGenomes, myTools


arguments = myTools.checkArgs(
    [("genesFile", file), ("transcriptsCoords", file)],
    [("useShortestTranscript", bool, True), ("sortOn5", bool, True), ("authorizedBiotypes", str, "protein_coding")],
    __doc__
    )

genome = myGenomes.Genome(arguments["genesFile"])
biotypes = set(arguments["authorizedBiotypes"].split(","))

# Loading transcripts list
lstTrans = collections.defaultdict(list)
f = myFile.myTSV.reader(arguments["transcriptsCoords"])
for l in f.csvobject:
    if l[-1] in biotypes:
        lstTrans[l[0]].append((int(l[2]), int(l[3]), l[1]))
f.file.close()

for chrom in genome.lstGenes:

    # list to sort
    tmp = []
    for gene in genome.lstGenes[chrom]:
        if arguments["useShortestTranscript"]:
            if gene.names[0] in lstTrans:
                l = [(y - x, x, y, t) for (x, y, t) in lstTrans[gene.names[0]]]
                best = min(l)
                if gene.strand == -1:
                    tmp.append((best[2], best[1]) + (gene, best[3]))
                else:
                    tmp.append(best[1:3] + (gene, best[3]))
            else:
                print("missing gene:", gene.names[0], file=sys.stderr)
        else:
            tmp.append((gene.beginning, gene.end, gene, None))

    #sort according to 5' or 3'
    if arguments["sortOn5"]:
        tmp.sort()
    else:
        import operator

        tmp.sort(key=operator.itemgetter(1))

    # stdout
    for (x, y, gene, name) in tmp:
        if gene.strand == -1:
            x, y = y, x
        res = [gene.chromosome, x, y, gene.strand, " ".join(gene.names)]
        if name is not None:
            res.append(name)
        print(myFile.myTSV.printLine(res))
