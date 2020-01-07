#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract the subtree descending from the first node validating the given condition"""
from LibsDyogen import myTools, myProteinTree


if __name__=='__main__':
    args = myTools.checkArgs([("forestFile", myTools.File), ("outFile", str),
                              ("id", int)], [], __doc__)

    node = args["id"]

    with open(args["outFile"], "w") as out:
        for tree in myProteinTree.loadTree(args["forestFile"]):
            if node in tree.info:
                tree.printTree(out, node)
                break
