#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Flatten, rebuild, and add family names."""

from __future__ import print_function
from sys import version_info, stdout, setrecursionlimit

if version_info[0] == 3:
    from LibsDyogen import myTools, myProteinTree, myPhylTree
else:
    from LibsDyogen.utils import myTools, myProteinTree, myPhylTree

from ToolsDyogen.treeTools.ALL.extractGeneFamilies import extractGeneFamilies


setrecursionlimit(10000)


def alwaysFalse(tree, node):
    return False


def alwaysTrue(tree, node):
    return True


def processTrees(ensemblTree, phylTree):
    for tree in myProteinTree.loadTree(arguments["ensemblTree"]):
        try:
            tree.flattenTree(phylTree, rec=True)
            # Not sure this step is useful, and why this hasLowScore function has no effect.
            tree.rebuildTree(phylTree, hasLowScore=alwaysFalse)
        except BaseException as err:
            err.args += ("Root id '%d'" % tree.root,)
            raise

        yield tree


if __name__ == '__main__':
    arguments = myTools.checkArgs(
        [("phylTree.conf", myTools.File), ("ensemblTree", myTools.File)],
        [("reuseNames", bool, False)],
        __doc__
    )


    phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
    setrecursionlimit(20000)

    count, dupCount, geneFamilies = extractGeneFamilies(phylTree,
                                            processTrees(arguments["ensemblTree"],
                                                         phylTree),
                                            arguments["reuseNames"])
