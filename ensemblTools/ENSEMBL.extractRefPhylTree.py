#!/usr/bin/env python3

"""
	Extrait l'arbre phylogenetique des especes utilise par Ensembl
"""

import sys
import io

from LibsDyogen import myFile, myTools, myPhylTree


arguments = myTools.checkArgs( [("IN.protein_tree_tag",file)], [], __doc__)

dicTaxonName = {}
dicTaxonID = {}
dicTaxonAlias = {}

# Chargement des donnees
print("Chargement des tags ...", end=' ', file=sys.stderr)
f = myFile.openFile(arguments["IN.protein_tree_tag"], "r")
for ligne in f:
	t = ligne[:-1].split("\t")
	if t[1] == "taxon_name":
		dicTaxonName[t[0]] = t[2]
	elif t[1] == "taxon_id":
		dicTaxonID[t[0]] = t[2]
	elif (t[1] == "taxon_alias") and (t[0] not in dicTaxonAlias):
		dicTaxonAlias[t[0]] = t[2]
	elif t[1] == "taxon_alias_mya":
		dicTaxonAlias[t[0]] = t[2]
	elif t[1] == "species_tree_string":
		tree = t[2]
print("OK (lengths:", len(dicTaxonName), len(dicTaxonID), len(dicTaxonAlias), ")", file=sys.stderr)

# Recoupement des infos
resTaxon = {}
for x in dicTaxonName:
	i = dicTaxonID[x]
	s = dicTaxonName[x]
	a = dicTaxonAlias.get(x, s)
	if x in resTaxon:
		# qui sont censees dire la meme chose
		assert resTaxon[i] == (s,a)
	else:
		resTaxon[i] = (s,a)
print(len(resTaxon), "taxa named", file=sys.stderr)

phylTree = myPhylTree.PhylogeneticTree(io.StringIO(tree))

# Impression sous mon format, avec des indentations
def do(node, indent):
	node = node.replace("*", "")
	print(("\t" * indent) + "|".join(resTaxon[node]))
	if node in phylTree.items:
		for (f,_) in phylTree.items[node]:
			do(f, indent+1)

do(phylTree.root, 0)

