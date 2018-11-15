#!/usr/bin/env python3

"""
	Blocs de syntenie entre deux especes
"""


import sys

from LibsDyogen import myPhylTree, myGenomes, myFile, myTools, myMaths, myDiags



# Arguments
modesOrthos = list(myDiags.OrthosFilterType._keys)
arguments = myTools.checkArgs( \
	[("genome1",file), ("genome2",file), ("ancGenes",file)], \
	[("fusionThreshold",int,-1), ("sameStrand",bool,True), ("orthosFilter",str,modesOrthos), ("minimalLength",int,2)], \
	__doc__ \
)

genome1 = myGenomes.Genome(arguments["genome1"])
genome2 = myGenomes.Genome(arguments["genome2"])
ancGenes = myGenomes.Genome(arguments["ancGenes"])
orthosFilter = myDiags.OrthosFilterType[modesOrthos.index(arguments["orthosFilter"])]

statsDiags = []
for ((c1,d1),(c2,d2),daa) in myDiags.calcDiags(genome1, genome2, ancGenes, \
	fusionThreshold=arguments["fusionThreshold"], sameStrand=arguments["sameStrand"], orthosFilter=orthosFilter, minChromLength=arguments["minimalLength"]):

	l = len(daa)
	if l < arguments["minimalLength"]:
		continue
	statsDiags.append(l)
	
	dic1 = genome1.lstGenes[c1]
	dic2 = genome2.lstGenes[c2]
	
	res = [l, \
		c1," ".join(genome1.lstGenes[c1][i1].names[0] for (i1,_) in d1), \
		c2," ".join(genome2.lstGenes[c2][i2].names[0] for (i2,_) in d2) ]
	
	res.append(myFile.myTSV.printLine(daa, " "))
	
	if arguments["sameStrand"]:
		# Un champ en plus pour l'orientation
		ds1 = [s1 for (_,s1) in d1]
		ds2 = [s2 for (_,s2) in d2]
		assert len(set(x/y for (x,y) in zip(ds1,ds2))) == 1, (ds1, ds2)
		res.append(myFile.myTSV.printLine(ds1, " "))

	print(myFile.myTSV.printLine(res))

print("Block length", myMaths.myStats.txtSummary(statsDiags), file=sys.stderr)

