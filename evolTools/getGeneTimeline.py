#! /usr/bin/env python

__doc__ = """
	Renvoie pour chaque gene ancestral le decompte des evenements qu'il subit sur chaque branche
"""

import sys

import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs(
    [("phylTree.conf", file), ("rootSpecies", str)],
    [("genesFiles", str, ""), ("ancGenesFiles", str, ""), ("countDup", bool, True), ("countLoss", bool, True)],
    __doc__
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Chargement des tous les fichiers
genes = {}
todo = {}
for e in phylTree.listSpecies:
    genes[e] = utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[e])
for a in phylTree.listAncestr:
    genes[a] = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[a])
    todo[a] = set(g.names[0] for g in genes[a])

allnames = set()
for a in phylTree.listAncestr:
    allnames.update(todo[a])
print(len(allnames), file=sys.stderr)

evol = dict((g, [None] * len(phylTree.indNames)) for g in allnames)


def countEvents(node, c, i):
    r = {}
    for (e, _) in phylTree.items.get(node, []):
        l = genes[e].getPositions(genes[node].lstGenes[c][i].names)
        if (len(l) > 1) and arguments["countDup"]:
            r[e] = 1
        elif (len(l) == 0) and (e not in phylTree.lstEsp2X) and arguments["countLoss"]:
            r[e] = 1
        else:
            r[e] = 0
        for (cc, ii) in l:
            r.update(countEvents(e, cc, ii))
    for (b, x) in r.items():
        evol[gene][phylTree.indNames[b]] = x
    return r


for gene in allnames:
    def lookupBeginning(node):
        if node in phylTree.listAncestr:
            if gene in todo[node]:
                (c, i) = genes[node].dicGenes[gene]
                countEvents(node, c, i)
                evol[gene][phylTree.indNames[node]] = 1
            else:
                for (e, _) in phylTree.items[node]:
                    lookupBeginning(e)


    lookupBeginning(arguments["rootSpecies"])
    print(gene, evol[gene])
