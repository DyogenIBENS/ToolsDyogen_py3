#!/usr/bin/env python3

"""
	Renvoie un arbre phylogenetique des especes avec les valeurs medianes issues des arbres de proteines
"""

import collections
import sys

from LibsDyogen import myTools, myPhylTree, myProteinTree


arguments = myTools.checkArgs([("phylTree.conf", file), ("proteinTree", file)], [], __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

lengths = collections.defaultdict(list)


# Parcours recursif de la famille de genes
def do(node):
    print("NEW TREE", file=sys.stderr)
    if node in tree.data:
        t1 = tree.info[node]['taxon_name']
        for (g, d) in tree.data[node]:
            # Une distance ne peut etre prise qu'entre deux noeuds de speciation
            if (tree.info[node]['Duplication'] == 0) and (tree.info[g]['Duplication'] == 0):
                t2 = tree.info[g]['taxon_name']
                # Les deux noeuds doivent etre strictement consecutifs
                if (phylTree.parent[t2].name == t1) and (d != 0):
                    lengths[(t1, t2)].append(d)
                    print(myFile.myTSV.printLine([t1, t2, d]), file=sys.stderr)
            do(g)


for tree in myProteinTree.loadTree(arguments["proteinTree"]):
    do(tree.root)

# On trie les listes des longueurs
for l in lengths.values():
    l.sort()


# Parcourt recursivement l'arbre et l'ecrit au format avec des parentheses, avec les longueurs de branche medianes
def convertToFlatFile(anc):
    a = phylTree.fileName[anc]
    if anc in phylTree.listSpecies:
        # On est arrive sur une feuille
        return a
    else:
        # On est sur un noeud, on construit la liste des distances
        l = []
        for (e, _) in phylTree.items[anc]:
            if len(lengths[(anc, e)]) == 0:
                # Par defaut, on revient a 0
                l.append(0)
            else:
                # La mediane est plus valable que la moyenne
                l.append(lengths[(anc, e)][len(lengths[(anc, e)]) / 2])
        # Constuction de la chaine finale
        return "(" + ",".join(
            [convertToFlatFile(e) + ":" + str(l) for ((e, _), l) in zip(phylTree.items[anc], l)]) + ")%s|%d" % (
        a, phylTree.ages[anc])


print(convertToFlatFile(phylTree.root), ";", file=sys.stdout)


def main():  # for setup.py console_scripts entry point
    pass
