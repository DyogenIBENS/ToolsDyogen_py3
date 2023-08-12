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
from LibsDyogen import myFile, myPhylTree, myProteinTree


def extractMultipleGeneTrees(proteinTree, family_name, field='family_name',
         toNewick=False, withAncSpeciesNames=False, withAncGenesNames=False,
         withTags=False, phyltree=None, output=None, force=False,
         mkdirs=False, firstmatch=False):
    if phyltree:
        phyltree = myPhylTree.PhylogeneticTree(phyltree)

    family_names = dict.fromkeys(family_name, 0)
    def get_tree_name(tree):
        return tree.info[tree.root].get('tree_name',
                tree.info[tree.root]['family_name'].split('.')[0])
    if field=='tree_name':
        get_field = get_tree_name
    else:
        def get_field(tree):
            return tree.info[tree.root][field]

    for tree in myProteinTree.loadTree(proteinTree):
        family = get_field(tree)
        if family in family_names:
            print("Found", family, end=' ', file=sys.stderr)
            wasfound = family_names[family]
            outfile = output.format(family=family, tree=get_tree_name(tree)) if output else '<stdout>'
            if os.path.isfile(outfile) and not wasfound and not firstmatch and not force:
                #if family_names[family] == 0:
                #FIXME so that you can omit the --force option but append to file
                print("%s exists. Skipping. (use --force)" % outfile, file=sys.stderr)
                family_names.pop(family)
            else:
                if phyltree is not None:
                    #markLowScore(tree, hasLowScore)
                    #flattenTree
                    #
                    tree.rebuildTree(phyltree)
                #TODO: start in new thread.
                filemode = 'a' if wasfound else 'w'
                try:
                    out = open(outfile, filemode) if output else sys.stdout
                except IOError:
                    if mkdirs:
                        os.makedirs(os.path.split(outfile)[0])
                        out = open(outfile, filemode)
                    else:
                        raise

                if toNewick:
                    print("Output to newick format", file=sys.stderr)
                    tree.printNewick(out, withDist=True, withTags=withTags,
                                     withAncSpeciesNames=withAncSpeciesNames,
                                     withAncGenesNames=withAncGenesNames,
                                     withID=withTags)
                else:
                    tree.printTree(out)
                if output: out.close()
                if firstmatch:
                    family_names.pop(family)
                else:
                    family_names[family] += 1
        if firstmatch and not family_names:
            break

    notfound = set((fam for fam,wasfound in family_names.items() if not wasfound))
    if notfound:
        print('WARNING: %d names were not found in field %r: %s' % (
              len(notfound), field, ' '.join(notfound)), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("proteinTree")
    parser.add_argument("family_name", nargs='+')
    parser.add_argument("-fromfile", action='store_true',
                        help='read `family_name` from a file, one per line.')
    parser.add_argument("-field", default="family_name",
                        choices=("tree_name", "family_name"),
                        help="[%(default)s]")
    parser.add_argument('-firstmatch', action='store_true',
                        help='Output only the first tree matching a family name.')
    parser.add_argument("-toNewick", action="store_true",
                        help="output in newick format")
    parser.add_argument("-withAncSpeciesNames", action="store_true")
    parser.add_argument("-noAncGenesNames", action="store_false",
                        dest='withAncGenesNames')
    parser.add_argument("-withTags", action="store_true")
    #parser.add_argument("-rebuild", action="store_true",
    #                    help="rebuild tree to fit species tree. Requires -phyltree")
    parser.add_argument("-phyltree", help=("path to PhylTree.conf file. -> rebuild"
                        "the gene tree to fit the species tree"))
    parser.add_argument("-output", "-o", #default='{genetree}.nwk',
                        help=("template for the filename. {tree} will be "
                              "replaced by the tree_name,"
                              " up to the first dot. {family} will be replaced"
                              " by the full family_name. [stdout]"))
    parser.add_argument("-force", "-f", action="store_true",
                        help="overwrite existing file")
    parser.add_argument("-mkdirs", action="store_true",
                        help="Create leading directories if needed")
    # TODO: add option to output protein_name in leaves, instead of gene_name
    # (need to add the option in LibsDyogen.MyProteinTree ...)

    arguments = parser.parse_args()
    if arguments.fromfile:
        fam_names = []
        for filename in arguments.family_name:
            with (sys.stdin if filename=='-' else open(filename)) as f:
                fam_names.extend(line.rstrip() for line in f if not line.startswith('#'))
        arguments.family_name = fam_names
    delattr(arguments, 'fromfile')

    extractMultipleGeneTrees(**vars(arguments))


if __name__=='__main__':
    main()
