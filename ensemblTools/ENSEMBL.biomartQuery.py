#! /usr/bin/env python


__doc__ = """
	Run the XMLfile BIOMART Query
	Usage:
		./ENSEMBL.biomartQuery.py XMLfiles/BIOMART.HumanProteinCodingGene.xml   -> will generate ouput.txt
		./ENSEMBL.biomartQuery.py XMLfiles/BIOMART.HumanProteinCodingGene.xml -outputFileName=HumanProteinCodingGene.txt
"""

import sys
import urllib.request, urllib.parse, urllib.error

import utils.myFile
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs( \
    [("xmlRequest", file)], \
    [("biomartServer", str, "http://www.ensembl.org/biomart/martservice"), ("outputFileName", str, "output.txt")], \
    __doc__ \
    )

# La requete
f = utils.myFile.openFile(arguments["xmlRequest"], "r")
request = f.read()
f.close()

print("Downloading XML Query", end=' ', file=sys.stderr)
urllib.request.urlretrieve(arguments["biomartServer"], filename=arguments["outputFileName"],
                   data=urllib.parse.urlencode({"query": request}))
print("OK", file=sys.stderr)
