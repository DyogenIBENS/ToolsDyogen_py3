#! /usr/bin/env python

"""
        From a species tree, print the number of extant species and ancetors.
	Optional: print the list of species
	
	Usage:	getInfoOnSpeciesTree.py PhylTree.conf
		getInfoOnSpeciesTree.py PhylTree.conf +speciesList +ancList
"""

import sys

from LibsDyogen import myFile, myTools, myPhylTree

arguments = myTools.checkArgs([("phylTree.conf", myTools.File)],
                              [("speciesList", bool, False),
                               ("ancList", bool, False)],
                              __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

if arguments["speciesList"]:
    print("Extant Species:", ",".join(x for x in phylTree.listSpecies), file=sys.stdout)
if arguments["ancList"]:
    print(file=sys.stdout)
    print("Ancestral Species:", ",".join(x for x in phylTree.listAncestr), file=sys.stdout)
print(file=sys.stdout)
print("Extant Species:", len(phylTree.listSpecies), file=sys.stdout)
print("Ancetral Species:", len(phylTree.listAncestr), file=sys.stdout)
