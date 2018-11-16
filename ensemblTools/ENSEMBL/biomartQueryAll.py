#!/usr/bin/env python3


"""
	Run the XMLQuery for each species in the PhylTree.conf file to the ENSEMBL biomart server
	Usage:
		./ENSEMBL.biomartQueryAll.py ../PhylTree.conf XMLfiles/BIOMART.transcriptsCoords.xml -outputFileName=coords.%s.list
		./ENSEMBL.biomartQueryAll.py ../PhylTree.conf XMLfiles/BIOMART.genesList.xml -outputFileName=genes.%s.list
"""

import sys
import time
import urllib.request, urllib.parse, urllib.error

from LibsDyogen import myFile, myPhylTree, myTools


# Arguments
arguments = myTools.checkArgs( \
    [("phylTree.conf", file), ("xmlRequest", file)], \
    [("biomartServer", str, "http://www.ensembl.org/biomart/martservice"), ("outputFileName", str, "output.%s.txt")], \
    __doc__ \
    )

# Read the Phylogenetic species tree
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Read the Query to execute
with myFile.openFile(arguments["xmlRequest"], "r") as f:
    request = f.read()

for esp in phylTree.listSpecies:
    # transform species name to ensembl species name  "Homo Sapiens" -> "hsapiens"
    tmp = esp.lower().split()
    tmp = tmp[0][0] + tmp[1]

    print("Downloading %s (%s) ..." % (esp, tmp), end=' ', file=sys.stderr)
    urllib.request.urlretrieve(arguments["biomartServer"], filename=arguments["outputFileName"] % phylTree.fileName[esp],
                       data=urllib.parse.urlencode({"query": request % tmp}))
    print("OK", file=sys.stderr)

    time.sleep(10)
