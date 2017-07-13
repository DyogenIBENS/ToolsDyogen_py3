#! /usr/bin/env python3

"""
Convert a species phylogenetic tree from/to newick or phyltree text format

usage:  ALL.convertNewickTree.py Species.conf -fromNewick > Species.nwk
        ALL.convertNewickTree.py Species.nwk +fromNewick > Species.conf
"""

import LibsDyogen.myFile     as myFile
import LibsDyogen.myTools    as myTools
import LibsDyogen.myPhylTree as myPhylTree

arguments = myTools.checkArgs([("phylTree.conf", myTools.File)],
                              [("fromNewick", bool, True)], __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

if arguments["fromNewick"]:

    # Returns the phyltree format (with indentation)
    def do(node, indent):
        node = node.replace("*", "")
        names = myFile.myTSV.printLine(
            [node] + [x for x in phylTree.commonNames.get(node, "")
                        if isinstance(x, str) and (x != node)], delim="|")
        print(("\t" * indent) + "%s" % names)
        if node in phylTree.items:
            for (f, _) in phylTree.items[node]:
                do(f, indent + 1)


    do(phylTree.root, 0)

else:
    # Returns the newick tree
    def convertToFlatFile(anc):

        a = phylTree.fileName[anc]  # anc.replace(' ', '.')
        if anc in phylTree.listSpecies:
            return a
        else:
            return "(" + ",".join(
                [convertToFlatFile(e) + ":" + str(l) for (e, l) in phylTree.items[anc]]) + ")%s|%d" % (
                a, phylTree.ages[anc])


    print(convertToFlatFile(phylTree.root), ";")
