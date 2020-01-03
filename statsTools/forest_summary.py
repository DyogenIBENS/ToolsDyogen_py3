#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Display some counts and summary statistics on a tree forest in the LibsDyogen format.
(ex: number of trees, of nodes, of leaves, of duplications, etc)

Usage:

    ./forest_summary.py [-h] <forestfile>

Options:
    -h  display this help.
"""

from sys import argv, stderr, exit
import numpy as np
import LibsDyogen.myProteinTree as ProteinTree


def forest_summary(forestfile):

    # Counts
    n_nodes = 0
    n_intnodes = 0
    n_leaves = 0
    n_dup = {}
    n_dubious = 0
    #n_multifurc = 0

    # Complete lists for summary stats
    taxa_set = set()
    species_set = set()
    dup_conf_scores = []
    dup_bootstraps = []
    
    tree_n_nodes = []
    tree_n_intnodes = []
    tree_n_leaves = []
    tree_n_speciesleaves = []

    tree_n_dup = []
    tree_n_dubious = []

    for tree_i, tree in enumerate(ProteinTree.loadTree(forestfile)):
        tree_n_intnodes.append(len(tree.data))
        tree_n_nodes.append(len(tree.info))
        
        # Counter for this tree
        tree_n_leaves.append(0)
        tree_n_speciesleaves.append(0)
        tree_n_dup.append(0)
        tree_n_dubious.append(0)

        for node_id, nodeinfo in tree.info.items():
            taxa_set.add(nodeinfo['taxon_name'])

            try:
                n_dup[nodeinfo['Duplication']] += 1
            except KeyError:
                n_dup[nodeinfo['Duplication']] = 1

            #if 'dubious_duplication' in nodeinfo:
            #    all_tree_n_dubious[-1] += 1
                
            if nodeinfo['Duplication'] != 0:
                tree_n_dup[-1] += 1
                dup_conf_scores.append(nodeinfo.get(
                                    'duplication_confidence_score',
                                    np.NaN))
                dup_bootstraps.append(nodeinfo.get(
                                    'Bootstrap',
                                    np.NaN))
                # dupli without conf scores: edited nodes with 'Duplication': 3
                # dupli without bootstrap: edited nodes with 'Duplication': 2

            if 'gene_name' in nodeinfo:
                species_set.add(nodeinfo['taxon_name'])
                tree_n_speciesleaves[-1] += 1

            if node_id not in tree.data:
                tree_n_leaves[-1] += 1

        assert tree_n_nodes[-1] - tree_n_intnodes[-1] == tree_n_leaves[-1]

    dup_conf_scores      = np.array(dup_conf_scores)
    dup_conf_scores_nan  = np.isnan(dup_conf_scores)
    dup_conf_scores      = dup_conf_scores[~dup_conf_scores_nan]
    dup_bootstraps       = np.array(dup_bootstraps)
    dup_bootstraps_nan   = np.isnan(dup_bootstraps)
    dup_bootstraps       = dup_bootstraps[~dup_bootstraps_nan]

    tree_n_nodes         = np.array(tree_n_nodes)
    tree_n_intnodes      = np.array(tree_n_intnodes)
    tree_n_leaves        = np.array(tree_n_leaves)
    tree_n_speciesleaves = np.array(tree_n_speciesleaves)

    tree_n_dup           = np.array(tree_n_dup)
    tree_n_dubious       = np.array(tree_n_dubious)

    return """
Nb of taxa    : {:d}
Nb of species : {:d}
Nb of trees   : {:d}

                       tot  tree average     tree std
n_nodes         : {:-8d}      {:-8.2f}     {:-8.2f}
n_intnodes      : {:-8d}      {:-8.2f}     {:-8.2f}
n_leaves        : {:-8d}      {:-8.2f}     {:-8.2f}
n_speciesleaves : {:-8d}      {:-8.2f}     {:-8.2f}

n_dup           : {}
                  tree average= {:4.2f}   tree std= {:4.2f}
#n_dubious
#n_multifurc

dup_conf_scores: average= {:8.5f}   std= {:8.5f}   missing= {:d}
dup_bootstraps:  average= {:8.5f}   std= {:8.5f}   missing= {:d}
""".format(
        len(taxa_set),
        len(species_set),
        (tree_i+1),
        tree_n_nodes.sum(),    tree_n_nodes.mean(),    tree_n_nodes.std(),
        tree_n_intnodes.sum(), tree_n_intnodes.mean(), tree_n_intnodes.std(),
        tree_n_leaves.sum(),   tree_n_leaves.mean(),   tree_n_leaves.std(),
        tree_n_speciesleaves.sum(), tree_n_speciesleaves.mean(), tree_n_speciesleaves.std(),
        ',  '.join('%d: %d' % item for item in n_dup.items()),
        tree_n_dup.mean(), tree_n_dup.std(),
        dup_conf_scores.mean(), dup_conf_scores.std(), dup_conf_scores_nan.sum(),
        dup_bootstraps.mean(), dup_bootstraps.std(), dup_bootstraps_nan.sum()
        )


def main():
    if len(argv) != 2:
        print('Bad number of arguments!\n' + __doc__, file=stderr)
        exit(1)
    elif argv[1] == '-h':
        print(__doc__, file=stderr)
        exit(0)
    print(forest_summary(argv[1]))


if __name__ == '__main__':
    main()
