#! /usr/bin/env python

__doc__ = """
	Print each pair of species and their commun ancestror
	
	Usage: misc.printLastCommonAncestors.py PhylTree.conf
"""

import sys
import itertools

import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs([("phylTree.conf", file)], [], __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

for a in phylTree.listAncestr:
    for (f1, f2) in itertools.combinations([f for (f, _) in phylTree.items[a]], 2):

        l1 = [e for e in phylTree.species[f1]]
        l2 = [e for e in phylTree.species[f2]]
        for (e1, e2) in itertools.product(l1, l2):
            print("%s\t%s\t%s" % (e1, e2, a), file=sys.stdout)
