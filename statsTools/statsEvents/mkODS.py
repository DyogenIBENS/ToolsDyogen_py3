#!/usr/bin/env python3

"""
	Cree le tableau des stats des blocs de syntenie et des genes ancestraux en fonction d'un seuil de coupure
"""

import sys

from LibsDyogen import myPhylTree, myGenomes, myFile, myTools, myMaths


def main():
    # Arguments
    arguments = myTools.checkArgs( \
        [("phylTree.conf", myTools.File), ("dirList", myTools.FileList(1))], \
        [("diagsFile", str, "diags/integr/diags.%s.list.bz2"), ("outputODS", str, "")], \
        __doc__ \
        )

    # L'arbre phylogenetique
    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


    # except KeyError:
    lstEspeces = sorted(set(phylTree.listAncestr))

    allCutoff = arguments["dirList"]

    titles = ["AncGenes", "Blocks", "Genes in blocks", "%Cov", "NbInt", "%CovInt", "Min", "25%", "50%", "75%", "N75", "N50",
              "N25", "Max", "Mean", "LongBlocks"]

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
            data[e].extend(myMaths.myStats.valSummary(lst)[:-2])

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
        from odfpy_examples import datatable

        textdoc = odf.opendocument.OpenDocumentSpreadsheet()

        for events in allEvents:
            # valevents = events.split("/")[-1]
            valevents = events
            # Premiere table avec les stats brutes
            val = [["Ancestor", "Age (My)"] + titles]
            for e in lstEspeces:
                val.append(alldata[events][e])

            table = datatable.DataTable(val)
            table.datasourcehaslabels = "both"
            t = table()
            t.setAttribute('name', valevents)
            textdoc.spreadsheet.addElement(t)


        # Table specifique pour un ancetre
        for esp in lstEspeces:
            # continue
            val = [["events"] + titles]
            for events in allEvents:
                # valevents = events.split("/")[-1]
                valevents = events
                val.append([valevents] + alldata[events][esp][2:])

            table = datatable.DataTable(val)
            table.datasourcehaslabels = "both"
            t = table()
            t.setAttribute('name', esp)
            textdoc.spreadsheet.addElement(t)


        # Resume final


        val = [["N50"] + ["events"] + [esp for esp in lstEspeces]]
        for events in allEvents:
            # valevents = events.split("/")[-1]
            valevents = events
            val.append([""] + [valevents] + [alldata[events][e][13] for e in lstEspeces])

        val.append(["Mean"] + ["events"] + [esp for esp in lstEspeces])
        for events in allEvents:
            # valevents = events.split("/")[-1]
            valevents = events
            val.append([""] + [valevents] + [alldata[events][e][16] for e in lstEspeces])

        val.append(["NbBlocks"] + ["events"] + [esp for esp in lstEspeces])
        for events in allEvents:
            # valevents = events.split("/")[-1]
            valevents = events
            val.append([""] + [valevents] + [alldata[events][e][3] for e in lstEspeces])
        val.append(["MaxLength"] + ["events"] + [esp for esp in lstEspeces])
        for events in allEvents:
            # valevents = events.split("/")[-1]
            valevents = events
            val.append([""] + [valevents] + [alldata[events][e][15] for e in lstEspeces])
        val.append(["LongBlocks"] + ["events"] + [esp for esp in lstEspeces])
        for events in allEvents:
            # valevents = events.split("/")[-1]
            valevents = events
            val.append([""] + [valevents] + [alldata[events][e][17] for e in lstEspeces])

        table = datatable.DataTable(val)
        table.datasourcehaslabels = "both"
        t = table()
        t.setAttribute('name', "Summary")
        textdoc.spreadsheet.addElement(t)

        textdoc.save(arguments["outputODS"])


if __name__ == '__main__':
    main()
