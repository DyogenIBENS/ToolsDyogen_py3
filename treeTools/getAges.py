#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Print the age of each taxon"""

from __future__ import print_function
from LibsDyogen import myTools, myPhylTree


if __name__ == '__main__':
    args = myTools.checkArgs([('phyltree', myTools.File)], [], __doc__)
    phyltree = myPhylTree.PhylogeneticTree(args['phyltree'])

    for taxon, age in sorted(phyltree.ages.items(), key=lambda x: (x[1], x[0])):
        print(taxon + '\t' + '%7g' % age)
