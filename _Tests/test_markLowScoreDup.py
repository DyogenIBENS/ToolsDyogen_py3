#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stderr, exit
from LibsDyogen import myProteinTree, myPhylTree
from treeTools.ENSEMBL.markLowScoreDup import countIndependentLosses


# Context
phyl_items = {'Hominoidea': [('Hominidae',1),  ('Nomascus leucogenys',1)],
              'Hominidae':  [('Homininae', 1), ('Pongo abelii', 1)],
              'Homininae':  [('HomoPan', 1),   ('Gorilla gorilla', 1)],
              'HomoPan':    [('Pan', 1),       ('Homo sapiens', 1)],
              'Pan':        [('Pan troglodytes', 1), ('Pan paniscus', 1)]}
phyl_officialNames = {name: name
                        for name in (set(phyl_items)
                                    | set(t for v in phyl_items.values()
                                            for t,_ in v))}

phyltree = myPhylTree.PhylogeneticTree((
                                        phyl_items,
                                        'Hominoidea',
                                        phyl_officialNames),
                                    skipInit=False,  # No effect if giving items
                                    stream=stderr)
phyltree.reinitTree(stream=stderr)


data = [
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),  (3, 0.1)]},
        info={1: {'Duplication': 2, 'taxon_name': 'Homo sapiens'},
              2: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              3: {'Duplication': 0, 'taxon_name': 'Homo sapiens'}},
        root=1),
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),  (3, 0.1)],
              2: [(4, 0.05), (5, 0.04)],
              3: [(6, 0.02), (7, 0.01)]},
        info={1: {'Duplication': 2, 'taxon_name': 'Pan'},
              2: {'Duplication': 0, 'taxon_name': 'Pan'},
              3: {'Duplication': 0, 'taxon_name': 'Pan'},
              4: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              5: {'Duplication': 0, 'taxon_name': 'Pan paniscus'},
              6: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              7: {'Duplication': 0, 'taxon_name': 'Pan paniscus'}},
        root=1),
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),  (3, 0.1)],
              2: [(4, 0.05), (5, 0.04)],
              3: [(6, 0.02), (7, 0.01)]},
        info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
              2: {'Duplication': 0, 'taxon_name': 'HomoPan'},
              3: {'Duplication': 0, 'taxon_name': 'HomoPan'},
              4: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              5: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              6: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              7: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'}},
        root=1),
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),  (3, 0.1)],
              2: [(4, 0.05), (5, 0.04)],
              3: [(6, 0.02), (7, 0.01)]},
        info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
              2: {'Duplication': 0, 'taxon_name': 'HomoPan'},
              3: {'Duplication': 0, 'taxon_name': 'HomoPan'},
              4: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              5: {'Duplication': 0, 'taxon_name': 'Pan paniscus'},
              6: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              7: {'Duplication': 0, 'taxon_name': 'Pan paniscus'}},
        root=1),
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),  (3, 0.1)]},
        info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
              2: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              3: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'}},
        root=1),
    #myProteinTree.ProteinTree(
    #    data={1: [(2, 0.1),  (3, 0.1)]},
    #    info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
    #          2: {'Duplication': 0, 'taxon_name': 'Pan paniscus'},
    #          3: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'}},
    #    root=1)
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),   (3, 0.1)],
              2: [(4, 0.05),  (5, 0.04)],
              3: [(6, 0.02),  (7, 0.01)],
              4: [(8, 0.005), (9, 0.005)],
              6: [(10, 0.005),(11, 0.005)]},
        info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
              2: {'Duplication': 0, 'taxon_name': 'HomoPan'},
              3: {'Duplication': 0, 'taxon_name': 'HomoPan'},
              4: {'Duplication': 0, 'taxon_name': 'Pan'},
              5: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              6: {'Duplication': 0, 'taxon_name': 'Pan'},
              7: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              8: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              9: {'Duplication': 0, 'taxon_name': 'Pan paniscus'},
              10: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              11: {'Duplication': 0, 'taxon_name': 'Pan paniscus'}},
        root=1),
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),   (3, 0.1)],
              3: [(6, 0.02),  (7, 0.01)],
              6: [(10, 0.005),(11, 0.005)]},
        info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
              2: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              3: {'Duplication': 0, 'taxon_name': 'HomoPan'},
              6: {'Duplication': 0, 'taxon_name': 'Pan'},
              7: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              8: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              9: {'Duplication': 0, 'taxon_name': 'Pan paniscus'},
              10: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              11: {'Duplication': 0, 'taxon_name': 'Pan paniscus'}},
        root=1),
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),   (6, 0.1)],
              6: [(10, 0.005),(11, 0.005)]},
        info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
              2: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              6: {'Duplication': 0, 'taxon_name': 'Pan'},
              10: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'},
              11: {'Duplication': 0, 'taxon_name': 'Pan paniscus'}},
        root=1),
    myProteinTree.ProteinTree(
        data={1: [(2, 0.1),   (10, 0.1)]},
        info={1: {'Duplication': 2, 'taxon_name': 'HomoPan'},
              2: {'Duplication': 0, 'taxon_name': 'Homo sapiens'},
              10: {'Duplication': 0, 'taxon_name': 'Pan troglodytes'}},
        root=1)
    ]

answers = [0, 0,   1./3, 1./2, 2./3, # 1.,]
               0, 0.5, 1, 1]

class Test_countIndependentLosses:

    def test_data0_1species_1dup(self):
        assert countIndependentLosses(data[0], data[0].root, phyltree) == answers[0]
    def test_data1_2species_1basaldup(self):
        assert countIndependentLosses(data[1], data[1].root, phyltree) == answers[1]
    def test_data2_3species_1basaldup_2del(self):
        assert countIndependentLosses(data[2], data[2].root, phyltree) == answers[2]
    def test_data3(self):
        assert countIndependentLosses(data[3], data[3].root, phyltree) == answers[3]
    def test_data4(self):
        assert countIndependentLosses(data[4], data[4].root, phyltree) == answers[4]
    def test_data5(self):
        assert countIndependentLosses(data[5], data[5].root, phyltree) == answers[5]
    def test_data6_3species_1basaldup_2postbasaldel(self):
        assert countIndependentLosses(data[6], data[6].root, phyltree) == answers[6]
    def test_data7_2species_1basaldup_2del(self):
        assert countIndependentLosses(data[7], data[7].root, phyltree) == answers[7]


def atestrunner_countIndependentLosses():
    # DEPRECATED. Use pytest.
    passed = []

    for i, (a, dat) in enumerate(zip(answers, data), start=1):
        result = countIndependentLosses(dat, dat.root, phyltree)

        if result == a:
            passed.append(1)
        else:
            passed.append(0)
            print('# %d. FAILED: %s (expected %s)' %(
                    i, result, a),
                  file=stderr)

    print('Passed %d/%d tests: [%s]' %(
            sum(passed), i, ' '.join(str(x) for x in passed)))

    return sum(passed) == i


if __name__ == '__main__':
    r = atestrunner_countIndependentLosses()
    exit(0 if r else 1)


