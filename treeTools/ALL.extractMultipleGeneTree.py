#! /usr/bin/env python3

"""
Search and save several gene trees from a forest of trees.

It can only search for "family_name" (faster as it's in the info of the root)
For "gene_name"/"protein_name", use `ALL.extractOneGeneTree.py`
"""

"""
USAGE:
    ./ALL.extractMultipleGeneTrees.py [options] <proteinTree> <family_name> [<family_name> ...]
"""

import sys
import os.path

import argparse
import LibsDyogen.myFile        as myFile
import LibsDyogen.myPhylTree    as myPhylTree
import LibsDyogen.myProteinTree as myProteinTree

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("proteinTree")
parser.add_argument("family_name", nargs='+')
#parser.add_argument("-field", default="gene_name",
#                    choices=("gene_name","family_name"))
parser.add_argument("-toNewick", action="store_true",
                    help="output in newick format")
parser.add_argument("-withAncSpeciesNames", action="store_true")
#parser.add_argument("-rebuild", action="store_true",
#                    help="rebuild tree to fit species tree. Requires -phyltree")
parser.add_argument("-phyltree", help=("path to PhylTree.conf file. -> rebuild"
                    "the gene tree to fit the species tree"))
parser.add_argument("-output", "-o",
                    help=("template for the filename. {genetree} will be "
                          "replaced by the name provided in the command line."))
parser.add_argument("-force", "-f", action="store_true",
                    help="overwrite existing file")
parser.add_argument("-mkdirs", action="store_true",
                    help="Create leading directories if needed")
# TODO: add option to output protein_name in leaves, instead of gene_name
# (need to add the option in LibsDyogen.MyProteinTree ...)

arguments = parser.parse_args()

if arguments.phyltree:
    phyltree = myPhylTree.PhylogeneticTree(arguments.phyltree)

family_names = set(arguments.family_name)

for tree in myProteinTree.loadTree(arguments.proteinTree):
    family = tree.info[tree.root]["family_name"]
    if family in family_names:
        print("Found", family, file=sys.stderr)
        outfile = arguments.output.format(genetree=family)
        if os.path.isfile(outfile) and not arguments.force:
            print("%s exists. Skipping. (use --force)" % outfile, file=sys.stderr)
        else:
            if arguments.phyltree:
                tree.rebuildTree(phyltree)
            try:
                out = open(outfile, 'w')
            except IOError:
                if arguments.mkdirs:
                    os.makedirs(os.path.split(outfile)[0])
                    out = open(outfile, 'w')
                else:
                    raise

            if arguments.toNewick:
                print("Output to newick format", file=sys.stderr)
                tree.printNewick(out, withDist=True, withTags=False,
                                 withAncSpeciesNames=arguments.withAncSpeciesNames,
                                 withAncGenesNames=True)
            else:
                tree.printTree(out)
            out.close()
        family_names.remove(family)
    if not family_names:
        break
