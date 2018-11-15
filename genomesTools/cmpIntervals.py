#!/usr/bin/env python3

"""
	Compute the conservation of each adjacency between extant and ancestral genome
	Usage:
		./cmpIntervals.py ../data/ancGenomes/genome.Boreoeutheria.list.bz2 ../data/genes/genesST.Homo.sapiens.list.bz2
"""

import sys
import itertools

from LibsDyogen import myTools, myGenomes



# Arguments:
arguments = myTools.checkArgs([("ancGenome", file), ("modernGenome", file)], [("minimalLength", int, 0)], __doc__)

ancGenome = myGenomes.Genome(arguments["ancGenome"])
genome = myGenomes.Genome(arguments["modernGenome"])


# Genome rewritting
def rewriteGenome(genome):
    newGenome = {}
    for chrom in genome.chrList[myGenomes.ContigType.Chromosome] + genome.chrList[
        myGenomes.ContigType.Scaffold]:
        if len(genome.lstGenes[chrom]) >= abs(arguments["minimalLength"]):
            newGenome[chrom] = [(gene.names[0], gene.strand) for gene in genome.lstGenes[chrom]]
    return newGenome


newGenome = rewriteGenome(genome)
print("modernGenome", sum([len(x) for x in newGenome.values()]), file=sys.stderr)

newAncGenome = rewriteGenome(ancGenome)
print("ancGenome", sum([len(x) for x in newAncGenome.values()]), file=sys.stderr)


# convertion tab between gene names (ancestral and modern)
def translate(genome1, genome1b, genome2, genome2b):
    trans = {}
    for gene in genome1:
        if gene.chromosome in genome1b:
            lpos = genome2.getPositions(gene.names)
            tmp = [genome2.lstGenes[pos.chromosome][pos.index].names[0] for pos in lpos if pos.chromosome in genome2b]
            if len(tmp) > 0:
                trans[gene.names[0]] = tmp
    return trans


translateMA = {}
for (x, y) in translate(genome, newGenome, ancGenome, newAncGenome).items():
    assert len(y) == 1
    translateMA[x] = y[0]
print("M->A", len(translateMA), list(translateMA.items())[0] if len(translateMA) > 0 else None, file=sys.stderr)

translateAM = translate(ancGenome, newAncGenome, genome, newGenome)
print("A->M", len(translateAM), list(translateAM.items())[0] if len(translateAM) > 0 else None, file=sys.stderr)

# assuming the two sets of genes are equivalents
assert set(translateMA.values()) == set(translateAM)
assert set(translateMA) == set(itertools.chain(*iter(translateAM.values())))


def removeNewSingletons(genome, translate):
    newGenome = {}
    for chrom in genome:
        tmp = [(g, s) for (g, s) in genome[chrom] if g in translate]
        if len(tmp) >= abs(arguments["minimalLength"]):
            newGenome[chrom] = genome[chrom]
    return newGenome


while arguments["minimalLength"] < -1:
    (n1, n2) = (len(newGenome), len(newAncGenome))
    print("iter", file=sys.stderr)
    newGenome = removeNewSingletons(newGenome, translateMA)
    print("modernGenome", sum([len(x) for x in newGenome.values()]), file=sys.stderr)
    newAncGenome = removeNewSingletons(newAncGenome, translateAM)
    print("ancGenome", sum([len(x) for x in newAncGenome.values()]), file=sys.stderr)

    if (n1, n2) == (len(newGenome), len(newAncGenome)):
        print("stop", file=sys.stderr)
        break

    translateMA = {}
    for (x, y) in translate(genome, newGenome, ancGenome, newAncGenome).items():
        assert len(y) == 1
        translateMA[x] = y[0]
    print("M->A", len(translateMA), list(translateMA.items())[0] if len(translateMA) > 0 else None, file=sys.stderr)

    translateAM = translate(ancGenome, newAncGenome, genome, newGenome)
    print("A->M", len(translateAM), list(translateAM.items())[0] if len(translateAM) > 0 else None, file=sys.stderr)

    # assuming the two sets of genes are equals
    assert set(translateMA.values()) == set(translateAM)
    assert set(translateMA) == set(itertools.chain(*iter(translateAM.values())))


# Intervals list
#     - all intervals
#     - Those with the two genes occurring in the other genome -> Eviction of specific genes
def listInterv(genome, translate):
    listIntAll = set()
    listIntFilt = []
    dicPos = {}
    dicLengths = {}
    for chrom in genome:
        listIntAll.update(myTools.myIterator.slidingTuple(genome[chrom]))
        tmp = [(g, s) for (g, s) in genome[chrom] if g in translate]
        for (i, (g, s)) in enumerate(tmp):
            dicPos[g] = (chrom, i, s)
        dicLengths[chrom] = len(tmp)
        listIntFilt.extend(myTools.myIterator.slidingTuple(tmp))
    print(len(listIntAll), len(listIntFilt), list(listIntAll)[0] if len(listIntAll) > 0 else None, \
    listIntFilt[0] if len(listIntFilt) > 0 else None, file=sys.stderr)
    return (listIntAll, listIntFilt, dicPos, dicLengths)


print("intMod", end=' ', file=sys.stderr)
(listIntMall, listIntMfilt, dicPosMod, dicModLengths) = listInterv(newGenome, translateMA)

print("intAnc", end=' ', file=sys.stderr)
(listIntAall, listIntAfilt, dicPosAnc, dicAncLengths) = listInterv(newAncGenome, translateAM)
listIntAfilt = set(listIntAfilt)
listIntAalls = set((g1, g2) for ((g1, s1), (g2, s2)) in listIntAall)
listIntAfilts = set((g1, g2) for ((g1, s1), (g2, s2)) in listIntAfilt)
assert len(listIntAalls) == len(listIntAall)
assert len(listIntAfilts) == len(listIntAfilt)

allendsT = ["NO_END", "ONE_END", "TWO_ENDS"]


def getRoom(length, pos, after):
    return (pos == (length - 1)) if after else (pos == 0)


# browsing from extant genome
seen = set()
for ((g1, s1), (g2, s2)) in listIntMfilt:
    tg1 = translateMA[g1]
    tg2 = translateMA[g2]
    flags = ["WITHOUT_GENE_GAIN_INSIDE" if ((g1, s1), (g2, s2)) in listIntMall else "WITH_GENE_GAIN_INSIDE"]

    if ((tg1, s1), (tg2, s2)) in listIntAfilt:
        status = "="
        flags.append("SAME_ORIENT")
        flags.append("WITHOUT_GENE_LOSS_INSIDE" if ((tg1, s1), (tg2, s2)) in listIntAall else "WITH_GENE_LOSS_INSIDE")
        seen.add((tg1, tg2))
    elif ((tg2, -s2), (tg1, -s1)) in listIntAfilt:
        status = "="
        flags.append("SAME_ORIENT")
        flags.append("WITHOUT_GENE_LOSS_INSIDE" if ((tg2, -s2), (tg1, -s1)) in listIntAall else "WITH_GENE_LOSS_INSIDE")
        seen.add((tg2, tg1))

    elif (tg1, tg2) in listIntAfilts:
        status = "="
        flags.append("DIFF_ORIENT")
        flags.append("WITHOUT_GENE_LOSS_INSIDE" if (tg1, tg2) in listIntAalls else "WITH_GENE_LOSS_INSIDE")
        seen.add((tg1, tg2))
    elif (tg2, tg1) in listIntAfilts:
        status = "="
        flags.append("DIFF_ORIENT")
        flags.append("WITHOUT_GENE_LOSS_INSIDE" if (tg2, tg1) in listIntAalls else "WITH_GENE_LOSS_INSIDE")
        seen.add((tg2, tg1))

    elif tg1 == tg2:
        status = "+"
        if s1 == s2:
            flags.append("DUPLICATES_SAME_ORIENT")
        else:
            flags.append("DUPLICATES_DIFF_ORIENT")

    else:
        status = "+"
        (ac1, ai1, as1) = dicPosAnc[tg1]
        (ac2, ai2, as2) = dicPosAnc[tg2]

        room1 = getRoom(dicAncLengths[ac1], ai1, as1 == s1)
        room2 = getRoom(dicAncLengths[ac2], ai2, as2 != s2)
        flags.append(allendsT[room1 + room2])

    print("\t".join([status, "%s/%d" % (g1, s1), "%s/%d" % (g2, s2), tg1, tg2] + flags))

# Browsing from ancestral genome (unconserved intervals)
for ((g1, s1), (g2, s2)) in listIntAfilt:
    if (g1, g2) in seen:
        continue
    status = "-"
    flags = ["WITHOUT_GENE_LOSS_INSIDE" if ((g1, s1), (g2, s2)) in listIntAall else "WITH_GENE_LOSS_INSIDE"]

    allends = []
    for (tg1, tg2) in itertools.product(translateAM[g1], translateAM[g2]):
        (ac1, ai1, as1) = dicPosMod[tg1]
        (ac2, ai2, as2) = dicPosMod[tg2]

        room1 = getRoom(dicModLengths[ac1], ai1, as1 == s1)
        room2 = getRoom(dicModLengths[ac2], ai2, as2 != s2)
        allends.append((room1 + room2, tg1, tg2))
    if len(allends) == 0:
        print(((g1, s1), (g2, s2)), translateAM[g1], translateAM[g2], file=sys.stderr)
    (i, tg1, tg2) = max(allends)
    flags.append(allendsT[i])

    print("\t".join([status, "%s/%d" % (g1, s1), "%s/%d" % (g2, s2), tg1, tg2] + flags))
