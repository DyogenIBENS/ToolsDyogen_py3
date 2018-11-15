#!/usr/bin/env python3

"""
	Renvoie les listes des devenirs de chaque gene le long des branches de l'arbre phylogenetique
"""

import sys

from LibsDyogen import myFile, myMaths, myTools, myGenomes, myPhylTree


# Arguments
arguments = myTools.checkArgs( [("phylTree.conf",file), ("genesFile",str), ("ancGenesFile",str)], [], __doc__)

# Chargement des tous les fichiers
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
genes = {}
for e in phylTree.listSpecies:
	genes[e] = myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
for a in phylTree.listAncestr:
	genes[a] = myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[a])

def transformName(esp, xxx_todo_changeme):
	(c,i) = xxx_todo_changeme
	if esp in phylTree.items:
		return i
	else:
		return str(c) + "|" + str(i)

def do(node):

	for (e,da) in phylTree.items.get(node, []):
		res = []
		seen = set([transformName(e,(c,i)) for (c,l) in genes[e].lstGenes.items() for i in range(len(l))])
		for g in genes[node].lstGenes[None]:
			lnewg = [transformName(e,x) for x in genes[e].getPositions(g.names)]
			seen.difference_update(lnewg)
			res.append(lnewg)
		#print (res,seen)
		nbPerdus = len([x for x in res if len(x) == 0])
		nbGagnes = len(seen)
		nbEgaux = len([x for x in res if len(x) == 1])
		nbFinal = sum([len(x) for x in res]) + nbGagnes
		nbDup = (nbFinal - nbGagnes) - (len(res) - nbPerdus)
		print(myFile.myTSV.printLine([node, e, da, len(res), nbFinal, nbPerdus, nbGagnes, nbEgaux, nbDup]))
		#print sum([len(x) for x in res])
		do(e)

#print (phylTree.root,len(genes[phylTree.root].lstGenes[None]))
print(myFile.myTSV.printLine(["parent", "fils", "Age", "nbInit", "nbFinal", "nbPerdus", "nbGagnes", "nbEgaux", "nbDup"]))
do(phylTree.root)

