#!/usr/bin/env python3

import sys
import collections

from LibsDyogen import myFile, myTools, myGenomes


arguments = myTools.checkArgs(
	[("genesFile",file), ("transcriptsCoords",file)],
	[("useShortestTranscript",bool,True), ("sortOn5",bool,True), ("authorizedBiotypes",str,"protein_coding")],
	"Cree une liste ordonnee des genes en tenant compte du plus petit transcrit"
)

genome = myGenomes.Genome(arguments["genesFile"])
biotypes = set(arguments["authorizedBiotypes"].split(","))

# Chargement de la liste des transcrits
lstTrans = collections.defaultdict(list)
f = myFile.myTSV.reader(arguments["transcriptsCoords"])
for l in f.csvobject:
	if l[-1] in biotypes:
		lstTrans[l[0]].append( (int(l[2]),int(l[3]),l[1]) )
f.file.close()

for chrom in genome.lstGenes:

	# Creation de la liste a trier
	tmp = []
	for gene in genome.lstGenes[chrom]:
		if arguments["useShortestTranscript"]:
			if gene.names[0] in lstTrans:
				l = [(y-x,x,y,t) for (x,y,t) in lstTrans[gene.names[0]]]
				best = min(l)
				tmp.append(best[1:3] + (gene,best[3]))
			else:
				print("missing gene:", gene.names[0], file=sys.stderr)
		else:
			tmp.append((gene.beginning,gene.end,gene,None))

	# Tri selon le 5' ou le 3'
	if arguments["sortOn5"]:
		tmp.sort()
	else:
		import operator
		tmp.sort(key=operator.itemgetter(1))

	# Affichage
	for (x,y,gene,name) in tmp:
		res = [gene.chromosome, x, y, gene.strand, " ".join(gene.names)]
		if name is not None:
			res.append(name)
		print(myFile.myTSV.printLine(res))

