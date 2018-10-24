#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""CLI Wrapper for all LibsDyogen.myProteinTree methods.

USAGE:

./forest.py <valid method> <proteinTreeFile> [<extra args>]

Example:

./forest.py flattenTree <proteinTreeFile> <PhylTreeFile> <recurs>
"""


from __future__ import print_function

from sys import argv, stderr, stdout, stdin
from CLItools.autoCLI import build_cli_processor

try:
    from LibsDyogen import myTools, myProteinTree, myPhylTree
except ImportError:
    try: 
        from LibsDyogen.utils import myTools, myProteinTree, myPhylTree
    except ImportError:
        from utils import myTools, myProteinTree, myPhylTree




def run(process, proteinTreeFile, converted_args):
    print("Give args: %s" % converted_args, file=stderr)

    if proteinTreeFile == '-':
        proteinTreeFile = stdin

    for tree in myProteinTree.loadTree(proteinTreeFile):
        process(tree, *converted_args)
        tree.printTree(stdout)
        


if __name__=='__main__':

    process, converted_args = build_cli_processor(myProteinTree.ProteinTree,
                                              {'phyltree': myPhylTree.PhylogeneticTree},
                                              1,
                                              *argv[1:])
    run(process, argv[2], converted_args)
