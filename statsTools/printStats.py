#! /usr/bin/env python


__doc__ = """
	Reads file of numbers and print statistics:
	Min  [Q25/Q50/Q75]  [N75/N50/N25]   Max   [Mean/Stddev-Length]
	
	
	Usage: ./printStats.py filename 
		./printStats.py filename +long +colNames
	
"""

import utils.myFile
import utils.myMaths
import utils.myTools

arguments = utils.myTools.checkArgs([("file", file)], [("long", bool, False), ("colNames", bool, False)], __doc__)

lst = []
f = utils.myFile.openFile(arguments["file"], 'r')

for l in f:
    c = l.split()
    for x in c:
        try:
            x = int(x)
        except ValueError:
            x = float(x)
        lst.append(x)
f.close()

# returns results

if arguments["long"]:
    if arguments["colNames"]:
        print(" ".join(("%s" % x) for x in
                       ["Min", "Q25", "Q50", "Q75", "N75", "N50", "N25", "WeightedAverage", "Max", "Mean", "Stddev",
                        "Length"]))

    print(" ".join(("%s" % x) for x in utils.myMaths.myStats.valSummary2(lst)))


else:

    print(utils.myMaths.myStats.syntheticTxtSummary(lst))
