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
                protein tree is rebuilt to match the species tree (UNTESTED!).
    -withAncSpeciesNames: if True, prepend the ancient species name to the
                          node name
"""

import sys

import LibsDyogen.myFile        as myFile
import LibsDyogen.myTools       as myTools
import LibsDyogen.myPhylTree    as myPhylTree
import LibsDyogen.myProteinTree as myProteinTree


def def_filtertest(fieldfilter, field='gene_name'):
    """Define the correct filter"""
    if field in ("gene_name", "protein_name", "tree_name"):
        def filtertest(nodeinfo):
            return nodeinfo.get(field) == fieldfilter
    elif field == "family_name":
        def filtertest(nodeinfo):
            return nodeinfo["family_name"].startswith(fieldfilter)
    else:
        raise RuntimeError("Invalid '-field' option")
    
    return filtertest


def search(filtertest, proteinTree, phyltree=None, toNewick=False,
           withAncSpeciesNames=False, withAncGenesNames=True):
    """search for the good gene tree, according to the function `filtertest`"""

    def do(node):
        if node in tree.data:
            for (g, d) in tree.data[node]:
                if do(g):
                    return True
        elif filtertest(tree.info[node]):
            return True
        return False

    for tree in myProteinTree.loadTree(proteinTree):
        if do(tree.root):
            if phyltree is not None:
                phyltree = myPhylTree.PhylogeneticTree(arguments['phyltree'])
                tree.rebuildTree(phyltree)
            if toNewick:
                print("Output to newick format", file=sys.stderr)
                tree.printNewick(sys.stdout, withDist=True, withTags=False,
                                 withAncSpeciesNames=withAncSpeciesNames,
                                 withAncGenesNames=withAncGenesNames)
            else:
                tree.printTree(sys.stdout)
            break


def main(**arguments):
    fieldfilter = arguments.pop("filter")
    field = arguments.pop("field")
    filtertest = def_filtertest(fieldfilter, field)
    search(filtertest, **arguments)


if __name__=='__main__':

    arguments = myTools.checkArgs([("proteinTree", myTools.File), ("filter", str)],
                                  [("phyltree", str, None),
                                   ("field", str, "gene_name"),
                                   ("toNewick", bool, False),
                                   ("withAncSpeciesNames", bool, False),
                                   ("withAncGenesNames", bool, True)],
                                  __doc__)
    main(**arguments)
