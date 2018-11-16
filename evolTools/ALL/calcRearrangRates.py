#!/usr/bin/env python3

"""
	Renvoie un tableau de statistiques de rearrangement le long des branches et un arbre de especes tenant compte de ces rearrangements.
"""

import sys

from LibsDyogen import myFile, myMaths, myTools, myGenomes, myPhylTree


# Argument:
arguments = myTools.checkArgs(
    [("phylTree.conf", file)],
    [("onlyOrthos", bool, False), ("in:genesFiles", str, ""), ("in:ancGenesFiles", str, ""), ("in:diagsFiles", str, ""),
     ("out:treeFile", str, "out.nwk"), ("out:statFile", str, "out.txt"), ("colNames", bool, True )],
    __doc__
)

# Chargement des tous les fichiers
###################################
phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

genes = {}
diags = {}
dicDiags = {}

for e in phylTree.listSpecies:
    # Les genes des especes modernes
    genes[e] = myGenomes.Genome(arguments["in:genesFiles"] % phylTree.fileName[e])
    diags[e] = []
    for (c, l) in genes[e].lstGenes.items():
        diags[e].append([((c, i), l[i].strand) for i in range(len(l))])

for a in phylTree.listAncestr:
    # Les genes ancestraux
    genes[a] = myGenomes.Genome(arguments["in:ancGenesFiles"] % phylTree.fileName[a])
    # Les diagonales
    diags[a] = []
    # On en profite pour lister les diagonales et les genes seuls
    notseen = set(range(len(genes[a].lstGenes[None])))
    f = myFile.openFile(arguments["in:diagsFiles"] % phylTree.fileName[a], "r")
    for l in f:
        t = l.split("\t")
        d = [int(x) for x in t[2].split()]
        s = [int(x) for x in t[3].split()]
        s = [2 * int(x >= 0) - 1 for x in s]
        diags[a].append(list(zip([(None, i) for i in d], s)))
        notseen.difference_update(d)
    f.close()
    assert len(notseen) == 0
# print >> sys.stderr, len(notseen)
# diags[a].extend( [((None,g),0)] for g in notseen)

# Creation des dictionnaires genes -> diags
for (esp, lst) in diags.items():
    dic = {}
    for (i, d) in enumerate(lst):
        for (j, (g, s)) in enumerate(d):
            dic[g] = (i, j, s)
    dicDiags[esp] = dic


# Les tables de conversion de diagonales entre ancetres successifs
val = {}


def do(node):
    # Les branches descendantes

    for (e, a) in phylTree.items.get(node, []):

        print("Stating branch %s -> %s ..." % (node, e), end=' ', file=sys.stderr)

        # Les correspondances de genes entre les deux noeuds
        corresp = {}
        for (c, i) in dicDiags[e]:
            l = genes[node].getPositions(genes[e].lstGenes[c][i].names)
            assert len(l) in [0, 1]
            if len(l) > 0:
                corresp[(c, i)] = dicDiags[node][tuple(l.pop())]
        print('e:', len(dicDiags[e]), 'node:', len(dicDiags[node]), 'corresp:', len(corresp), file=sys.stderr)
        nbOK = 0
        nbNO = 0
        nbXX = 0
        for d in diags[e]:
            if arguments["onlyOrthos"]:
                d = [(g, s) for (g, s) in d if g in corresp]
            for ((g1, s1), (g2, s2)) in myTools.myIterator.slidingTuple(d):
                if (g1 not in corresp) or (g2 not in corresp):
                    continue
                (c1, i1, t1) = corresp[g1]
                (c2, i2, t2) = corresp[g2]
                j1 = i2 - s2 * t2
                j2 = i1 + s1 * t1
                print("XXX:", (s1, s2), corresp[g1], corresp[g2], (
                    j1, j2), i2 == j2, i1 == j1, s1 * s2 == t1 * t2, len(
                    diags[node][c1]), len(diags[node][c2]), end=' ', file=sys.stderr)
                if c1 == c2:
                    if (i2 == j2) and (s1 * s2 == t1 * t2):
                        nbOK += 1
                        print("OK1", file=sys.stderr)
                    else:
                        nbNO += 1
                        print("NO1", file=sys.stderr)
                elif (j2 >= 0) and (j2 < len(diags[node][c1])):
                    print("NO2", file=sys.stderr)
                    nbNO += 1
                else:
                    j1 = i2 - s2 * t2
                    if (j1 >= 0) and (j1 < len(diags[node][c2])):
                        nbNO += 1
                        print("NO3", file=sys.stderr)
                    else:
                        nbXX += 1
                        print("??", file=sys.stderr)

        print(myFile.myTSV.printLine(
            [node, e, nbOK, nbNO, nbXX, nbOK + nbNO + nbXX, (100. * nbOK) / (nbOK + nbNO), len(diags[node]), a]), file=outfileTxt)
        print("OK", file=sys.stderr)
        val[(node, e)] = float(nbNO) / (nbOK + nbNO)
        do(e)


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
            l.append(val[(anc, e)])
        # Constuction de la chaine finale
        return "(" + ",".join(
            [convertToFlatFile(e) + ":" + str(l) for ((e, age), l) in zip(phylTree.items[anc], l)]) + ")%s|%d" % (
            a, phylTree.ages[anc])


# ecriture du fichier de stats
outfileTxt = myFile.openFile(arguments["out:statFile"], "w")

if (arguments("colNames")):
    print(myFile.myTSV.printLine(
        ["anc", "desc", "nbOK", "nbNO", "nbXX", "nbOK+nbNO", "nbNO+nbXX", "(100. * nbOK) / (nbOK + nbNO)", "nbDiags_Anc",
        "DivTime MA"]), file=outfileTxt)
do(phylTree.root)
outfileTxt.close()


# ecriture de l'arbre.
outfileTree = myFile.openFile(arguments["out:treeFile"], "w")
print(convertToFlatFile(phylTree.root), ";", file=outfileTree)
outfileTree.close()
