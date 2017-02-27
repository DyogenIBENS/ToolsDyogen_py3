#!/usr/bin/python

__doc__ = """
	Find Strong ancGenes Families 1-1 (no duplication, no loss in descendants.
"""

import sys

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("target",str), ("IN.ancGenesFiles",str), ("OUT.ancGenesFiles",str)], [("except2XSpecies",bool,True)], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
target = phylTree.officialName[arguments["target"]]

if arguments["except2XSpecies"] == "True":

    lstAncGenomes=[x for x in phylTree.listAncestr if phylTree.dicParents[x][target] == target and x not in phylTree.lstEsp2X]
    lstModernGenomes=[x for x in phylTree.listSpecies if phylTree.dicParents[x][target] == target and x not in phylTree.lstEsp2X]

else:
    lstAncGenomes=[x for x in phylTree.listAncestr if phylTree.dicParents[x][target] == target]
    lstModernGenomes=[x for x in phylTree.listSpecies if phylTree.dicParents[x][target] == target]


print("loading ancestral genomes", lstAncGenomes, file=sys.stderr)
phylTree.loadSpeciesFromList([x for x in lstAncGenomes], arguments["IN.ancGenesFiles"])
print("loading Modern genomes", lstModernGenomes, file=sys.stderr)
phylTree.loadSpeciesFromList([x for x in lstModernGenomes], arguments["IN.ancGenesFiles"])

#loading gene->genome dictionnary

extantGenes={}

for x in lstModernGenomes:
    for (i,gene) in enumerate(phylTree.dicGenomes[x].lstGenes[None]):
        extantGenes[gene.names[1]]=x

#print >> sys.stderr, extantGenes

for x in lstAncGenomes:
        if arguments["except2XSpecies"]=="True":
            lstDescSpecies = [y for y in phylTree.listSpecies if phylTree.dicParents[y][x] == x and y not in phylTree.lstEsp2X]
        else:
            lstDescSpecies = [y for y in phylTree.listSpecies if phylTree.dicParents[y][x] == x]

        if len(lstDescSpecies) > 0:
            f= utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[x], "w")

            print(x, lstDescSpecies, file=sys.stdout)
            for ancGene in phylTree.dicGenomes[x]:
                #print >> sys.stderr, ancGene
                nbDesc={}
                ancGenename=ancGene.names[0]
                for descSpecies in lstDescSpecies:
                   nbDesc[descSpecies] = 0
                for modernGene in ancGene.names[1:]:

                       if (modernGene in extantGenes):
                            #print >> sys.stderr, "modernGene:", modernGene,extantGenes[modernGene]
                            nbDesc[extantGenes[modernGene]]+=1
                       else:
                           next
                #print >> sys.stderr, nbDesc
                a=list(nbDesc.values())
                #print a, max(a), min(a)
                if  max(a)>1 or min(a)==0:
                    print(ancGenename, file=f)
                else:
                    print(" ".join(x for x in ancGene.names ), file=f)

            f.close()