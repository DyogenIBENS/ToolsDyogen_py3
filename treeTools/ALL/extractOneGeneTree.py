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


def do(node, filtertest):
    """filtertest: takes `node.info` in input, return True/False if matches
    your search criterion."""
    if node in tree.data:
        for (g, d) in tree.data[node]:
            if do(g):
                return True
    elif filtertest(tree.info[node]):
        return True
    return False


def search(filtertest, proteinTree, phyltree=None, toNewick=False, withAncSpecieNames=False):
    """search for the good gene tree, according to the function `filtertest`"""
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


if __name__=='__main__':

    arguments = myTools.checkArgs([("proteinTree", myTools.File), ("filter", str)],
                                  [("phyltree", str, None),
                                   ("field", str, "gene_name"),
                                   ("toNewick", bool, False),
                                   ("withAncSpeciesNames", bool, False)],
                                  __doc__)

    field = arguments.pop("field")
    fieldfilter = arguments.pop("filter")

    # Define the correct filter
    if field in ("gene_name", "protein_name"):
        def filtertest(nodeinfo):
            return nodeinfo[field] == fieldfilter
    elif field == "family_name":
        def filtertest(nodeinfo):
            return nodeinfo["family_name"].startswith(fieldfilter)
    else:
        print("Invalid '-field' option", file=sys.stderr)
        sys.exit(1)


    search(filtertest, **arguments)
