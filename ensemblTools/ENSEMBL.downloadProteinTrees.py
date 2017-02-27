#! /usr/bin/env python

__doc__ = """
	Met en forme les arbres de proteines (en les telechargeant / utilisant les fichiers locaux)
	Format gene trees in phyltree format with downloaded sql files from Ensembl. (old version, use the perl script instead)
"""

import os
import sys
import collections

import utils.myFile
import utils.myTools
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-XXX/mysql/ensembl_compara_XXX"), \
	("IN.member",str,"member.txt.gz"), \
	("IN.genome_db",str,"genome_db.txt.gz"), \
	("IN.protein_tree_node",str,"protein_tree_node.txt.gz"), \
	("IN.protein_tree_member",str,"protein_tree_member.txt.gz"), \
	("IN.protein_tree_tag",str,"protein_tree_tag.txt.gz"), \
	("IN.protein_tree_stable_id",str,"protein_tree_stable_id.txt.gz"), \
	("OUT.tree",str,"tree.1.ensembl"), \
	("OUT.geneInfo",str,"geneInfoFromTrees.txt"), \
	], \
	__doc__ \
)


phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


##############################################
# Chargement de la base de donnees d'Ensembl #
##############################################


# On charge les liens taxon_id -> species name
###############################################
print("Chargement des liens taxon_id -> species_name ...", end=' ', file=sys.stderr)
taxonName = {}
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.genome_db"]), "r")
for ligne in utils.myFile.myTSV.MySQLFileLoader(f):
	t = ligne.split("\t")
	taxonName[t[1]] = t[2]
f.close()
print(len(taxonName), "especes OK", file=sys.stderr)


# On charge les liens member_id -> protein name
################################################
print("Chargement des liens member_id -> protein_name ...", end=' ', file=sys.stderr)
tmpLinks = {}
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.member"]), "r")
for ligne in utils.myFile.myTSV.MySQLFileLoader(f):
	t = ligne.split("\t")
	# A un numero member_id, on associe les noms (gene/transcrit/proteine) et l'espece
	if t[7] != "\\N":
		x = t[8].split()
		tmpLinks[t[0]] = ((x[1].split(':')[1], x[0].split(':')[1], t[1]), taxonName[t[4]])
f.close()
del taxonName
print(len(tmpLinks), "membres OK", file=sys.stderr)


# On charge les liens node_id -> member_id
###########################################
print("Chargement des liens node_id -> member_id ...", end=' ', file=sys.stderr)
x = 0
info = collections.defaultdict(dict)
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_member"]), "r")
fi = utils.myFile.openFile(arguments["OUT.geneInfo"], "w")
for ligne in utils.myFile.myTSV.MySQLFileLoader(f):
	t = ligne.split("\t")
	data = tmpLinks[t[2]]
	info[int(t[0])] = {'gene_name': data[0][0], 'taxon_name': data[1]}
	print(utils.myFile.myTSV.printLine([t[0], data[1], data[0][0], data[0][1], data[0][2], t[4]]), file=fi)
	x += 1
f.close()
fi.close()
print(x, "proteines OK", file=sys.stderr)
del tmpLinks


# On charge les liens node_id -> tree_name
###########################################
print("Chargement des node_id -> tree_name ...", end=' ', file=sys.stderr)
x = 0
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_stable_id"]), "r")
for ligne in utils.myFile.myTSV.MySQLFileLoader(f):
	t = ligne.split("\t")
	info[int(t[0])]["tree_name"] = t[1]
	x += 1
f.close()
print(x, "noms OK", file=sys.stderr)


# On charge les liens node_id -> infos 
#######################################
print("Chargement des liens node_id -> infos ...", end=' ', file=sys.stderr)
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_tag"]), "r")
for ligne in utils.myFile.myTSV.MySQLFileLoader(f):
	t = ligne.split("\t")

	# On extrait l'info sous forme numerique si possible
	try:
		t[2] = int(t[2])
	except ValueError:
		try:
			t[2] = float(t[2])
		except ValueError:
			pass
	
	info[ int(t[0]) ][t[1]] = t[2]

f.close()
print(len(info), "infos OK", file=sys.stderr)


# On charge les liens node_id_father -> node_id_son
####################################################
print("Chargement des arbres (node_id_father -> node_id_son) ...", end=' ', file=sys.stderr)
data = collections.defaultdict(list)
f = utils.myFile.openFile(os.path.join(arguments["IN.EnsemblURL"], arguments["IN.protein_tree_node"]), "r")
for ligne in utils.myFile.myTSV.MySQLFileLoader(f):
	t = ligne.split("\t")
	# On rajoute le fils a son pere (avec une distance)
	node = int(t[1])
	data[ node ].append( (int(t[0]), float(t[6])) )
f.close()
allroots = set(root for (root,_) in data[1])
print(len(data), "branches OK", file=sys.stderr)


# On regle les parametres de info
##################################
acceptedTagNames = set(["Bootstrap", "Duplication", "duplication_confidence_score", "dubious_duplication", "taxon_name", "gene_name"])
for (node,inf) in info.items():
	if 'taxon_name' in inf:
		# Pour passer des '/' et ' ' a '-' et '_', et enlever le 'Silurana'
		inf['taxon_name'] = phylTree.officialName[inf['taxon_name']]
	if node not in allroots:
		for tagname in list(inf.keys()):
			if tagname not in acceptedTagNames:
				del inf[tagname]

for dat in data.values():
	dat.sort()

# On a besoin des genomes modernes pour reconnaitre les genes
print("Mise en forme des arbres ...", end=' ', file=sys.stderr)
ft = utils.myFile.openFile(arguments["OUT.tree"], "w")
for root in allroots:
	if 'taxon_name' in info[root]:
		utils.myProteinTree.ProteinTree(data, info, root).printTree(ft)
ft.close()
print(len(data[1]), "arbres OK", file=sys.stderr)

