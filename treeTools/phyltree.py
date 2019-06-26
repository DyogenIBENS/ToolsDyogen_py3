#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Print PhylTree attribute values."""


from __future__ import print_function

import argparse as ap
from LibsDyogen import myPhylTree


def main(phyltreefile, attr=None, keys=None):
    phyltree = myPhylTree.PhylogeneticTree(phyltreefile)
    if not attr:
        output = ('Available attributes:\n' +
                  '\n'.join('%-25s %s' %(a, type(getattr(phyltree,a)))
                            for a in sorted(phyltree.__dict__)
                            if not a.startswith('_') and not callable(getattr(phyltree,a))))
        print(output)
        return
    
    value = getattr(phyltree, attr)
    if isinstance(value, (str, int, float)):
        print(value)
        return

    if keys:
        def getvalue(key, default=''):
            try:
                return value[key]
            except KeyError:
                return ''
        try:
            values = [(k, getvalue(k)) for k in keys]
        except TypeError:
            def getvalue(key):
                return key in value  # what if value is a string?
            values = [(k, getvalue(k)) for k in keys]

        if len(keys) > 1:
            output = '\n'.join('%s\t%s' % item for item in values)
            print(output)
            return
        else:
            value = values[0][1]

    try:
        output = '\n'.join('%s\t%s' % (k,v) for k,v in sorted(value.items()))
    except (AttributeError,TypeError):
        try:
            output = '\n'.join(str(v) for v in sorted(value))
        except TypeError:
            output = '%s\n%s' % (type(value), value)

    print(output)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('phyltreefile')
    parser.add_argument('attr', nargs='?')
    parser.add_argument('keys', nargs='*')
    
    args = parser.parse_args()
    main(**vars(args))

