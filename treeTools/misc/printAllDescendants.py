#! /usr/bin/env python

"""
	Print all the ancestral and/or extants descendants species from one Ancestor

	Usage:
	    misc.printAllDescendants.py PhylTree.conf Boreoeutheria +withExtantSpecies
	    misc.printAllDescendants.py PhylTree.conf Percomorpha
"""

import sys

from LibsDyogen import myTools, myPhylTree

arguments = myTools.checkArgs([("phylTree.conf", myTools.File),
                               ("taxon_name", str)],
                              [("withExtantSpecies", bool, True)], __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

fromAnc = arguments["taxon_name"]

descendantAnc=[]
extantSpecies=[]
for desc in phylTree.allDescendants[fromAnc]:
    if (len(phylTree.allDescendants[desc]) > 1):
        descendantAnc.append(desc)
    else:
        extantSpecies.append(desc)


if ( arguments["withExtantSpecies"]):

    for anc in descendantAnc:
        print(anc, file=sys.stdout)

    print("----------------------", file=sys.stdout)
    for esp in extantSpecies:
        print(esp, file=sys.stdout)

else:
    for anc in descendantAnc:
        print(anc, file=sys.stdout)
