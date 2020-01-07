#!/usr/bin/env python3



from LibsDyogen import myTools, myPhylTree, myProteinTree


def pruneProteinTree(tree, node, subtree):
    if tree.info[node]['taxon_name'] not in subtree


def main(forestfile, phyltreefile, speciesfile, invert=False):
    phyltree = myPhylTree.PhylogeneticTree(phyltreefile)
    with open(speciesfile) as f:
        badspecies = [line.rstrip() for line in f if not line.startswith('#')]

    subroot, subtree = phyltree.getSubTree(badspecies)

    for tree in myProteinTree.loadTree(forestfile):



if __name__ == '__main__':

    args = myTools.checkArgs([("phyltreefile", myTools.File),
                              ("forestfile",   myTools.File),
                              ("speciesfile",   myTools.File)],
                             [("invert", bool, False)],
                             __doc__)

    main(**args)
