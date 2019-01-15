#!/usr/bin/env python3


"""
	Run the XMLfile BIOMART Query
	Usage:
		./ENSEMBL.biomartQuery.py XMLfiles/BIOMART.HumanProteinCodingGene.xml   -> will generate ouput.txt
		./ENSEMBL.biomartQuery.py XMLfiles/BIOMART.HumanProteinCodingGene.xml -outputFileName=HumanProteinCodingGene.txt
"""

from __future__ import print_function

import sys
import urllib.request, urllib.parse, urllib.error

from LibsDyogen import myFile, myTools

# Arguments
arguments = myTools.checkArgs(
    [("xmlRequest", myTools.File)],
    [("biomartServer", str, "http://www.ensembl.org/biomart/martservice"),
     ("outputFileName", str, "output.txt")],
    __doc__)

# La requete
with myFile.openFile(arguments["xmlRequest"], "r") as f:
    request = f.read()

print("Downloading XML Query", end=' ', file=sys.stderr)
urllib.request.urlretrieve(arguments["biomartServer"], filename=arguments["outputFileName"],
                   data=urllib.parse.urlencode({"query": request}).encode())
print("OK", file=sys.stderr)
