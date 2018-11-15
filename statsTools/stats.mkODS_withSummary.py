#!/usr/bin/env python3

"""
	
	Creates a spreadsheet ods of stats for each ancestral reconstruction according to threshold
	Gives a summary for each threshold:
	
	NbAncGenes/N50/WA/Moyenne/NbBlocks/L70

	Compute gains.

	Usage: 	stats.mkODS_withSummary.py PhylTree.Ensembl.82.conf trees/?.??/diags/integr/denovo-size-custom.refine-all.extend-all.halfinsert-all/ -diagsFile=anc/diags.%s.list.bz2 -outputODS=stats_nogroups_1.ods
		stats.mkODS_withSummary.py PhylTree.Ensembl.82.conf trees/?.??/diags/integr/final/ -diagsFile=anc/diags.%s.list.bz2 -outputODS=stats_final.ods

"""

import sys

from LibsDyogen import myPhylTree, myGenomes, myFile, myTools, myMaths


# Arguments
arguments = myTools.checkArgs( \
    [("phylTree.conf", file), ("dirList", myTools.FileList(1))], \
    [("diagsFile", str, "diags/integr/final/anc/diags.%s.list.bz2"), ("outputODS", str, "")], \
    __doc__ \
    )

# L'arbre phylogenetique
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
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
          "N25", "WeigthedAverage", "Max", "Mean", "LongBlocks"]

alldata = {}
alldiff = {}
allEvents = []

for cutoff in allCutoff:
    # allEvents.append(cutoff.replace(".refine32-all.fuseSingletons-all.halfInsert-all.groups","").replace("denovo-",""))
    allEvents.append(cutoff)
for events in allEvents:

    print(events, "...", end=' ', file=sys.stderr)

    # Recuperation des donnees de longueur de blocs
    alldata[events] = data = {}
    for e in lstEspeces:
        # print >> sys.stderr, e, "...",
        f = myFile.openFile(events + "/" + (arguments["diagsFile"] % phylTree.fileName[e]), "r")
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
        data[e].extend(myMaths.myStats.valSummary2(lst)[:-2])

        # on trie la liste des blocks par taille de blocks.
        lstSort = list(lst)
        lstSort.sort()
        # print  >> sys.stderr, lst
        nbBlock = 0
        ValKaryo75 = (tot - sing) * 75 / 100
        Karyo75 = 0
        while Karyo75 < ValKaryo75:
            tmp = lstSort.pop()
            Karyo75 += tmp
            nbBlock += 1

        data[e].append(nbBlock)
        print(e, "...", nbBlock, "...", end=' ', file=sys.stderr)
    if events == allEvents[0]:
        ref = data
    # else:

    # alldiff[events] = diff = {}
    #	for e in lstEspeces:
    #		newdata = [(x-ref[e][i] if i >= 2 else x) for (i,x) in enumerate(data[e])]
    #		newdata.insert(2, 100*(1.-float(newdata[4])/newdata[2]) if newdata[2] != 0 else None)
    #		diff[e] = newdata
    print("OK", file=sys.stderr)

if arguments["outputODS"] == "":
    for events in allEvents:
        print(events, file=sys.stdout)
        print(myFile.myTSV.printLine(["Ancestor", "Age (My)"] + titles))
        for e in lstEspeces:
            print(myFile.myTSV.printLine(alldata[events][e]))
    if events in alldiff:
        print(myFile.myTSV.printLine(["Ancestor", "Age (My)", "%Useful Gene Loss"] + titles))
        for e in lstEspeces:
            print(myFile.myTSV.printLine(alldiff[events][e]))

else:
    import odf.opendocument
    import odfpy.datatable

    textdoc = odf.opendocument.OpenDocumentSpreadsheet()

    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        # Premiere table avec les stats brutes
        val = [["Ancestor", "Age (My)"] + titles]
        for e in lstEspeces:
            val.append(alldata[events][e])

        table = odfpy.datatable.DataTable(val)
        table.datasourcehaslabels = "both"
        t = table()
        t.setAttribute('name', valevents)
        textdoc.spreadsheet.addElement(t)

    # if events in alldiff:
    #
    #			# Deuxieme table avec les differences par rapport a la reference
    #			val = [["Ancestor", "Age (My)", "%Useful Gene Loss"] + titles]
    #			for e in lstEspeces:
    #				val.append(alldiff[events][e])
    #
    #			table = odfpy.datatable.DataTable(val)
    #			table.datasourcehaslabels = "both"
    #			t = table()
    #			t.setAttribute('name', "d"+valevents)
    #			textdoc.spreadsheet.addElement(t)

    # Table specifique pour un ancetre
    for esp in lstEspeces:
        # continue
        val = [["events"] + titles]
        for events in allEvents:
            # valevents = events.split("/")[-1]
            valevents = events
            val.append([valevents] + alldata[events][esp][2:])

        table = odfpy.datatable.DataTable(val)
        table.datasourcehaslabels = "both"
        t = table()
        t.setAttribute('name', esp)
        textdoc.spreadsheet.addElement(t)


    # Resume final
    # val = [["events", "Mean gain", "Median gain", "N50 gain", "%Cov gain", "%CovInt gain", "BlockLength %gain (mean)", "BlockLength %gain (Median)", "BlockLength %gain (N50)", "Cov %gain", "CovInt %gain"]]
    #	for events in allEvents:
    #		valevents = events.split("/")[-1]
    #		val.append( [valevents] + [myMaths.myStats.mean([alldiff[events][e][i] for e in lstEspeces]) for i in [17, 12, 14, 6, 8]] +
    #			[myMaths.myStats.mean([100*float(alldata[events][e][i-1]-alldata[allEvents[0]][e][i-1])/alldata[allEvents[0]][e][i-1] for e in lstEspeces]) for i in [17, 12, 14, 6, 8]]
    #		)
    #	table = odfpy.datatable.DataTable(val)
    #	table.datasourcehaslabels = "both"
    #	t = table()
    #	t.setAttribute('name', "events")
    #	textdoc.spreadsheet.addElement(t)

    # Pour les courbes
    val = [["AncGenes"] + ["events"] + [esp for esp in lstEspeces]]
    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        val.append([""] + [valevents] + [alldata[events][e][2] for e in lstEspeces])

    val.append(["WeigthedAverage"] + ["events"] + [esp for esp in lstEspeces])
    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        val.append([""] + [valevents] + [int(alldata[events][e][15]) for e in lstEspeces])

    val.append(["N50"] + ["events"] + [esp for esp in lstEspeces])

    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        val.append([""] + [valevents] + [alldata[events][e][13] for e in lstEspeces])

    val.append(["Mean"] + ["events"] + [esp for esp in lstEspeces])
    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        val.append([""] + [valevents] + [alldata[events][e][17] for e in lstEspeces])

    val.append(["NbBlocks"] + ["events"] + [esp for esp in lstEspeces])
    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        val.append([""] + [valevents] + [alldata[events][e][3] for e in lstEspeces])
    val.append(["MaxLength"] + ["events"] + [esp for esp in lstEspeces])
    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        val.append([""] + [valevents] + [alldata[events][e][16] for e in lstEspeces])
    val.append(["LongBlocks"] + ["events"] + [esp for esp in lstEspeces])
    for events in allEvents:
        # valevents = events.split("/")[-1]
        valevents = events
        val.append([""] + [valevents] + [alldata[events][e][18] for e in lstEspeces])

    table = odfpy.datatable.DataTable(val)
    table.datasourcehaslabels = "both"
    t = table()
    t.setAttribute('name', "Summary")
    textdoc.spreadsheet.addElement(t)

    textdoc.save(arguments["outputODS"])
