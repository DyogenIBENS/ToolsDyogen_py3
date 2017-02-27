#! /usr/bin/env python3

"""
Search and print out a single gene tree from the forest of trees.

USAGE:
    ./ALL.extractOneGeneTree.py GeneTrees.bz2 ENSOARG00000011785

OPTIONS:
    -field    : field on which to select: possible values:
                    - gene_name     : modern gene ID (default)
                    - protein_name  : modern protein ID
                    - family_name : use a gene tree ID instead : 'ENSGT...'
    -toNewick : output in newick format
    -phyltree : a PhylTree.conf file (species phylogeny). If given, the 
                protein tree is rebuilt to match the species tree.
    -withAncSpeciesNames: if True, prepend the ancient species name to the
                          node name
"""

import sys

import LibsDyogen.myFile        as myFile
import LibsDyogen.myTools       as myTools
import LibsDyogen.myPhylTree    as myPhylTree
import LibsDyogen.myProteinTree as myProteinTree

arguments = myTools.checkArgs([("proteinTree", myTools.File), ("filter", str)],
                              [("phyltree", str, None),
                               ("field", str, "gene_name"),
                               ("toNewick", bool, False),
                               ("withAncSpeciesNames", bool, False)],
                              __doc__)

field = arguments["field"]

if field in ("gene_name", "protein_name"):
    def filtertest(nodeinfo):
        return nodeinfo[field] == arguments["filter"]
elif field == "family_name":
    def filtertest(nodeinfo):
        return nodeinfo["family_name"].startswith(arguments["filter"])
else:
    print("Invalid '-field' option", file=sys.stderr)
    sys.exit(1)


def do(node):
    if node in tree.data:
        for (g, d) in tree.data[node]:
            if do(g):
                return True
    elif filtertest(tree.info[node]):
        return True
    return False


# searching for the good gene tree
for tree in myProteinTree.loadTree(arguments["proteinTree"]):
    if do(tree.root):
        if arguments['phyltree']:
            phyltree = myPhylTree.PhylogeneticTree(arguments['phyltree'])
            tree.rebuildTree(phyltree)
        if arguments['toNewick']:
            print("Output to newick format", file=sys.stderr)
            tree.printNewick(sys.stdout, withDist=True, withTags=False,
                             withAncSpeciesNames=arguments['withAncSpeciesNames'],
                             withAncGenesNames=True)
        else:
            tree.printTree(sys.stdout)
        break
