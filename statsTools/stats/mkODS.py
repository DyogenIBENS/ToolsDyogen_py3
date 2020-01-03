#!/usr/bin/env python3

"""

	Creates a spreadsheet ods of stats for each ancestral reconstruction according to threshold
	
	Compute gains.

	Usage: 	stats.mkODS.py PhylTree.Ensembl.82.conf trees/?.??/diags/integr/denovo-size-custom.refine-all.extend-all.halfinsert-all/ -diagsFile=anc/diags.%s.list.bz2 -outputODS=stats_nogroups_1.ods
		stats.mkODS.py PhylTree.Ensembl.82.conf trees/?.??/diags/integr/final/ -diagsFile=anc/diags.%s.list.bz2 -outputODS=stats_final.ods
"""

import sys

from LibsDyogen import myPhylTree, myGenomes, myFile, myTools, myMaths


def main():
    # Arguments
    arguments = myTools.checkArgs( \
        [("phylTree.conf", myTools.File), ("dirList", myTools.FileList(1))], \
        [("diagsFile", str, "diags/integr/final/anc/diags.%s.list.bz2"), ("outputODS", str, "")], \
        __doc__ \
        )

    # L'arbre phylogenetique
    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

    # Liste des especes dans le bon ordre
    todo = set(phylTree.listAncestr)
    try:
        l1 = phylTree.dicLinks["Euteleostomi"]["Homo sapiens"][:-1]
        todo.difference_update(l1)
        l2 = phylTree.dicLinks["Glires"]["Murinae"]
        todo.difference_update(l2)
        l3 = [e for e in todo if phylTree.isChildOf(e, "Mammalia")]
        l3 = sorted(l3, key=lambda e: phylTree.ages[e], reverse=True)
        todo.difference_update(l3)
        l4 = [e for e in todo if phylTree.isChildOf(e, "Clupeocephala")]
        l4 = sorted(l4, key=lambda e: phylTree.ages[e], reverse=True)
        todo.difference_update(l4)
        l5 = [e for e in todo if phylTree.isChildOf(e, "Amniota")]
        l5 = sorted(l5, key=lambda e: phylTree.ages[e], reverse=True)
        todo.difference_update(l5)
        l6 = sorted(todo, key=lambda e: phylTree.ages[e], reverse=True)
        lstEspeces = l6 + l5 + l4 + l1 + l3 + l2
    except KeyError:
        lstEspeces = sorted(phylTree.listAncestr)
    # lstEspeces = l5

    # lstEspeces = ["Euteleostomi", "Amniota", "Boreoeutheria"]

    allCutoff = arguments["dirList"]

    titles = ["AncGenes", "Blocks", "Genes in blocks", "%Cov", "NbInt", "%CovInt", "Min", "25%", "50%", "75%", "N75", "N50",
              "N25", "Max", "Mean"]

    alldata = {}
    alldiff = {}

    for cutoff in allCutoff:

        print(cutoff, "...", end=' ', file=sys.stderr)

        # Recuperation des donnees de longueur de blocs
        alldata[cutoff] = data = {}
        for e in lstEspeces:

            f = myFile.openFile(cutoff + "/" + (arguments["diagsFile"] % phylTree.fileName[e]), "r")
            lst = []
            sing = 0
            tot = 0
            interv = 0
            for l in f:
                x = int(l.split("\t")[1])
                tot += x
                if x >= 2:
                    lst.append(x)
                    interv += (x - 1)
                else:
                    sing += 1
            f.close()

            data[e] = [e, phylTree.ages[e], tot, len(lst), tot - sing, (100. * (tot - sing)) / tot, interv,
                       (100. * interv) / (tot - 20.)]
            data[e].extend(myMaths.myStats.valSummary(lst)[:-2])

        if cutoff == allCutoff[0]:
            ref = data
        # else:

        alldiff[cutoff] = diff = {}
        for e in lstEspeces:
            newdata = [(x - ref[e][i] if i >= 2 else x) for (i, x) in enumerate(data[e])]
            newdata.insert(2, 100 * (1. - float(newdata[4]) / newdata[2]) if newdata[2] != 0 else None)
            diff[e] = newdata
        print("OK", file=sys.stderr)

    if arguments["outputODS"] == "":
        for cutoff in allCutoff:
            print(myFile.myTSV.printLine(["Ancestor", "Age (My)"] + titles))
            for e in lstEspeces:
                print(myFile.myTSV.printLine(alldata[cutoff][e]))
        if cutoff in alldiff:
            print(myFile.myTSV.printLine(["Ancestor", "Age (My)", "%Useful Gene Loss"] + titles))
            for e in lstEspeces:
                print(myFile.myTSV.printLine(alldiff[cutoff][e]))

    else:
        import odf.opendocument
        from odfpy_examples import datatable

        textdoc = odf.opendocument.OpenDocumentSpreadsheet()

        for cutoff in allCutoff:
            valCutoff = cutoff.split("/")[-1]

            # Premiere table avec les stats brutes
            val = [["Ancestor", "Age (My)"] + titles]
            for e in lstEspeces:
                val.append(alldata[cutoff][e])

            table = datatable.DataTable(val)
            table.datasourcehaslabels = "both"
            t = table()
            t.setAttribute('name', valCutoff)
            textdoc.spreadsheet.addElement(t)

            if cutoff in alldiff:

                # Deuxieme table avec les differences par rapport a la reference
                val = [["Ancestor", "Age (My)", "%Useful Gene Loss"] + titles]
                for e in lstEspeces:
                    val.append(alldiff[cutoff][e])

                table = datatable.DataTable(val)
                table.datasourcehaslabels = "both"
                t = table()
                t.setAttribute('name', "d" + valCutoff)
                textdoc.spreadsheet.addElement(t)

        # Table specifique pour un ancetre
        for esp in lstEspeces:
            # continue
            val = [["cutoff"] + titles]
            for cutoff in allCutoff:
                valCutoff = cutoff.split("/")[-1]
                val.append([valCutoff] + alldata[cutoff][esp][2:])

            table = datatable.DataTable(val)
            table.datasourcehaslabels = "both"
            t = table()
            t.setAttribute('name', esp)
            textdoc.spreadsheet.addElement(t)


        # Resume final
        val = [["cutoff", "Mean gain", "Median gain", "N50 gain", "%Cov gain", "%CovInt gain", "BlockLength %gain (mean)",
                "BlockLength %gain (Median)", "BlockLength %gain (N50)", "Cov %gain", "CovInt %gain"]]
        for cutoff in allCutoff:
            valCutoff = cutoff.split("/")[-1]
            val.append([valCutoff] + [myMaths.myStats.mean([alldiff[cutoff][e][i] for e in lstEspeces]) for i in
                                      [17, 12, 14, 6, 8]] +
                       [myMaths.myStats.mean([100 * float(
                           alldata[cutoff][e][i - 1] - alldata[allCutoff[0]][e][i - 1]) / alldata[allCutoff[0]][e][i - 1]
                                                    for e in lstEspeces]) for i in [17, 12, 14, 6, 8]]
                       )
        table = datatable.DataTable(val)
        table.datasourcehaslabels = "both"
        t = table()
        t.setAttribute('name', "cutoff")
        textdoc.spreadsheet.addElement(t)

        textdoc.save(arguments["outputODS"])
