#!/usr/bin/env python3

"""
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys
import collections

from LibsDyogen import myDiags, myMaths, myTools, myGenomes, myPhylTree




# Argument:
arguments = myTools.checkArgs( \
    [("phylTree.conf", file)], \
    [("IN.genesFile", str, ""), ("IN.ancGenesFile", str, ""), ("IN.diagsFile", str, "")], \
    __doc__ \
    )

# Chargement des tous les fichiers
###################################
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

genes = {}
diags = {}
dicDiags = {}

for e in phylTree.listSpecies:
    # Les genes des especes modernes
    genes[e] = myGenomes.Genome(arguments["IN.genesFile"] % phylTree.fileName[e])
    diags[e] = [[g] for g in range(len(list(genes[e])))]

for a in phylTree.listAncestr:
    # Les genes ancestraux
    genes[a] = myGenomes.Genome(arguments["IN.ancGenesFile"] % phylTree.fileName[a])
    # Les diagonales
    tmp = myDiags.loadDiagsFile(arguments["IN.diagsFile"] % phylTree.fileName[a], [a], phylTree.officialName)[a]
    # On en profite pour lister les diagonales et les genes seuls
    notseen = set(range(len(genes[a].lstGenes[None])))
    diags[a] = []
    for (d, _, _, _, _) in tmp:
        notseen.difference_update(d)
        diags[a].append(d)
    diags[a].extend([g] for g in notseen)

# Creation des dictionnaires genes -> diags
for (esp, lst) in diags.items():
    dic = {}
    for (i, d) in enumerate(lst):
        for (j, g) in enumerate(d):
            dic[g] = (i, j)
    dicDiags[esp] = dic
    print((esp, lst))


# Les tables de conversion de diagonales entre ancetres successifs
def do(node):
    # Les branches descendantes
    for (e, _) in phylTree.items.get(node, []):

        print("Stating branch %s -> %s ..." % (node, e), end=' ', file=sys.stderr)

        # Les correspondances de genes entre les deux noeuds
        corresp = [[j for (_, j) in genes[node].getPosition(g.names)] for g in genes[e]]

        res = []
        for d in diags[e]:
            count = collections.defaultdict(int)
            for g in d:
                for gg in corresp[g]:
                    if gg in dicDiags[node]:
                        count[dicDiags[node][gg][0]] += 1
            res.append(list(count.items()))
        print(res)
        print("OK", file=sys.stderr)
        do(e)


do(phylTree.root)
