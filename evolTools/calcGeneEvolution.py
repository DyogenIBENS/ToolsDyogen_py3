#! /usr/bin/env python

__doc__ = """
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs([("phylTree.conf", file), ("rootSpecies", str)],
                                    [("genesFile", str, ""), ("ancGenesFile", str, "")], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
# Chargement des tous les fichiers
genes = {}
for e in phylTree.listSpecies:
    genes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
for a in phylTree.listAncestr:
    genes[a] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[a])


def transformName(esp, xxx_todo_changeme):
    (c, i) = xxx_todo_changeme
    return genes[esp].lstGenes[c][i].names[0]


def do(node):
    for (e, _) in phylTree.items.get(node, []):
        res = {}
        seen = set([transformName(e, (c, i)) for (c, l) in genes[e].lstGenes.items() for i in range(len(l))])
        for g in genes[node].lstGenes[None]:
            lnewg = [transformName(e, x) for x in genes[e].getPositions(g.names)]
            seen.difference_update(lnewg)
            res[g.names[0]] = lnewg
        print((res, seen))
        do(e)


print([gene.names[0] for gene in genes[arguments["rootSpecies"]]])

do(arguments["rootSpecies"])
