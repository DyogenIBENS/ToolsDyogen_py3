#! /usr/bin/env python

__doc__ = """
	Extrait (des genomes reels) la liste des evenements de duplications/pertes/gains sur chaque branche de l'arbre
"""

import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree

arguments = utils.myTools.checkArgs([("phylTree.conf", file)],
                                    [("rootSpecies", str, ""), ("genesFile", str, ""), ("ancGenesFile", str, "")],
                                    __doc__)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


@utils.myTools.memoize
def getGenome(e):
    if e in phylTree.listSpecies:
        return utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
    else:
        return utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[e])


def transformName(esp, xxx_todo_changeme):
    (c, i) = xxx_todo_changeme
    return getGenome(esp).lstGenes[c][i].names[0]


def do(node):
    for (e, _) in phylTree.items.get(node, []):
        trans = {}
        seen = set([transformName(e, (c, i)) for (c, l) in getGenome(e).lstGenes.items() for i in range(len(l))])
        deleted = set()
        for g in getGenome(node).lstGenes[None]:
            lnewg = [transformName(e, x) for x in getGenome(e).getPositions(g.names)]
            seen.difference_update(lnewg)
            if len(lnewg) > 0:
                trans[g.names[0]] = lnewg
            else:
                deleted.add(g.names[0])

        print(utils.myFile.myTSV.printLine([node, e, len(trans), len(deleted), len(seen), trans, deleted, seen]))
        #print utils.myFile.myTSV.printLine([node, e, deleted])
        do(e)


root = phylTree.root if len(arguments["rootSpecies"]) == 0 else arguments["rootSpecies"]
print(utils.myFile.myTSV.printLine([root, set(gene.names[0] for gene in getGenome(root))]))
do(root)
