#!/usr/bin/env python3

"""
	Parcourt un fichier de genome et enleve les genes inclus dans un autre
"""

import sys
import collections
import itertools
import operator

from LibsDyogen import myFile, myTools, myGenomes


# Arguments
arguments = myTools.checkArgs( [("genome",file)], [], __doc__)

genome = myGenomes.Genome(arguments["genome"])

for c in genome.lstGenes:

	lref = list(genome.lstGenes[c])
	lref.sort(key=operator.attrgetter("beginning"))

	lnew = list(genome.lstGenes[c])
	lnew.sort(key=operator.attrgetter("end"))

	comb = myTools.myCombinator()
	for (g1,g2) in zip(lref, lnew):
		if g1 != g2:
			comb.addLink([g1,g2])

	removed = set()
	for grp in comb:
		grp.sort(key=operator.attrgetter("beginning"))
		while len(grp) > 0:
			first = grp.pop(0)
			toremove = [x for x in grp if x.end <= first.end]
			removed.update(toremove)
			grp = [x for x in grp if x not in removed]

	for gene in lref:
		if gene not in removed:
			print(myFile.myTSV.printLine([c, gene.beginning, gene.end, gene.strand, " ".join(gene.names)]))

