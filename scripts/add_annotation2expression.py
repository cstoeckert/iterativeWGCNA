#!/usr/bin/env python
# pylint: disable=invalid-name
'''
Adds gene annotation to expression data file
prints result to stdout
'''

from __future__ import print_function

import sys
import argparse

def parse_annotation():
    '''
    parse annotation file and create annotation map
    '''
    with open(args.annotation, 'r') as f:
        f.next()
        for line in f:
            values = line.rstrip().split('\t')
            annotationMap[values[0]] = values[1]


def write_annotated_expression():
    '''
    write updated expression file with new annotation added in
    '''
    with open(args.expression, 'r') as f:
        header = f.next().rstrip().split('\t')

        # find index of mapping column
        mappingIndex = header.index(args.mappingColumn)

        # update header
        header.insert(1, args.name)

        print('\t'.join(header), file=sys.stdout)

        for line in f:
            values = line.rstrip().split('\t')
            key = values[mappingIndex]
            annotation = 'UNCLASSIFIED' if key == 'UNCLASSIFIED' else annotationMap[key]
            values.insert(1, annotation)
            print('\t'.join(values), file=sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='add_annotation2expression.py',
                                     description="adds annotation to gene expression file",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-e', '--expression',
                        metavar='<gene expression file>',
                        help="full path to input gene expression file; "
                        + "if full path is not provided,\n"
                        + "assumes the file is in the current directory\n;"
                        + "NOTE: assumes header for all files\n",
                        required=True)

    parser.add_argument('-a', '--annotation',
                        metavar='<annotation file>',
                        help="full path to file containing annotation\n"
                        + "expects 2 tab-delim columns (key\tvalue)\n"
                        + "NOTE: assumes file has header\n",
                        required=True)

    parser.add_argument('-m', '--mappingColumn',
                        metavar='<column name>',
                        help="column name used for mapping in the expression dataset",
                        required=True)

    parser.add_argument('-n', '--name',
                        metavar='<annotation name>',
                        help="column name for annotation in updated expression dataset",
                        required=True)

    args = parser.parse_args()
    annotationMap = {}

    try:
        parse_annotation()
        write_annotated_expression()
    except IOError as e:
        sys.exit("I/O error({0}): {1}".format(e.errno, e.strerror))
    except ValueError as e:
        sys.exit("ValueError: Field " + args.name + " not in expression dataset")

