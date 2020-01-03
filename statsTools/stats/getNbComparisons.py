#!/usr/bin/env python3

"""

	print the number of pairwise comparison for each ancestor. 

	Usage:	stats.getNbComparisons.py PhylTree.conf -diags=$GENOMICUS/data82/trees/0.30/diags/integr/final/anc/diags.%s.list.bz2 +colNames
		stats.getNbComparisons.py PhylTree.conf 

"""

import sys
import itertools

from LibsDyogen import myFile, myTools, myMaths, myPhylTree


def main():
    # Arguments
    arguments = myTools.checkArgs( \
        [("phylTree.conf", myTools.File)], [("diags", str, ""), ("colNames", bool, False)], \
        __doc__ \
        )

    # L'arbre phylogenetique
    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

    if (arguments["colNames"]):
        print(myFile.myTSV.printLine(
            ["Ancestor", "NbComp", "Nb(In/Out)Comp", "Nb(In/In)Comp", "Age", "MeanSize_OfBlocks","N50Size_OfBlocks", "WASize_OfBlocks", "NbComp/Age"]), file=sys.stdout)

    for anc in phylTree.listAncestr:
        # nb d'outgroup:
        ###############

        nb_outgroup = len(phylTree.outgroupSpecies[anc])

        # nb d'Ingroups.
        ##############
        nbInSpec = [len(phylTree.species[x]) for (x, _) in phylTree.items[anc]]

        l = [len(phylTree.species[x]) for (x, _) in phylTree.items[anc]]
        # for (x,_) in phylTree.items[anc]:
        #	print >> sys.stderr, phylTree.species[x]
        l.append(nb_outgroup)

        # Comp InSpecies/OutGroups
        #########################

        compInOut = sum(nb_outgroup * n1 for n1 in nbInSpec)

        # Comp InSpecies/InSpecies
        #########################

        compInIn = sum(n1 * n2 for (n1, n2) in itertools.combinations(nbInSpec, 2))

        nbc = sum(n1 * n2 for (n1, n2) in itertools.combinations(l, 2))

        # quid des blocs.
        ###############

        totalStat= []
        if (arguments["diags"] != ""):
            r = []
            f = myFile.openFile(arguments["diags"] % phylTree.fileName[anc], "r")
            for line in f:
                x = int(line.split("\t")[1])
                if x > 1:
                    r.append(x)
            f.close()
            #lll = float(sum(r)) / len(r)
            totalStat = myMaths.myStats.valSummary2(r)
        else:
            lll = "NONE"
        ###############

        print(myFile.myTSV.printLine(
            [anc, nbc, compInOut, compInIn, phylTree.ages[anc], totalStat[9], totalStat[6], int(totalStat[7]), float(nbc) / phylTree.ages[anc]]))


if __name__ == '__main__':
    main()
