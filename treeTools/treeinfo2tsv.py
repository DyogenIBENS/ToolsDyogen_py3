#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import LibsDyogen.myProteinTree as ProteinTree


def main(forestfile, outputfile):
    header = ['id', 'root', 'family_name', 'Duplication', 'duplication_confidence_score', 'Bootstrap', ]
    with open(outputfile, 'r') as out:
