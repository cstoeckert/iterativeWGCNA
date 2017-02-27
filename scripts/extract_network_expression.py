#!/usr/bin/env python

'''
Given original expression data (iterativeWGCNA input)
and final membership file from iterativeWGCNA output,
extracts the expression data for a selected module

If the option inclMembership is specified, will insert a
column for module assignments between gene labels
and expression values
'''

from __future__ import print_function

import argparse
from sys import stdout

def parse_id2gene_map():
    '''
    returns a dict of numeric IDs to gene labels
    if geneIdMapping was provided
    '''
    geneIdMap = None
    if args.geneIdMapping is not None:
        geneIdMap = {}
        with open(args.geneIdMapping, 'r') as f:
            f.next() # skip header
            for line in f:
                values = line.strip().split('\t')
                geneIdMap[values[0]] = values[1]
 

def parse_membership():
    '''
    returns a dict of gene identifiers to module
    membership assignments
    '''
    membership = {}
    with open(args.classification, 'r') as f:
        f.next() # skip header
        for line in f:
            values = line.strip().split('\t')
            membership[values[0]] = values[1]


def extract_expression():
    '''
    extracts expression according to selected options
    and prints to stdout
    '''
    with open(args.expression, 'r') as f:
        fields = f.next().strip().split('\t')
        if args.addMembershipColumn:
            fields.insert(1, 'Module')
        print('\t'.join(header), file=stdout)

        for line in f:
            expression = line.strip().split('\t')
            geneId = expression.pop()
            geneMembership = membership[gene]
            
            if geneMembership == 'UNCLASSIFIED' \
                and args.includeUnclassified is None:
                continue
            
            if (args.module is not None \
              and geneMembership == args.module) \
              or args.module is None:

              # map id
              if args.addMembershipColumn:
                expression.insert(0, geneMembership)

                print('\t'.join(
              


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='extract_network_expression.py',
                                     description="extracts expression for specified module and performs mapping between expression data and module membership",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-e', '--expression',
                        metavar='<gene expression file>',
                        help="full path to input gene expression file; "
                        + "if full path is not provided,\n"
                        + "assumes the file is in the current directory\n;"
                        + "file format same as expected for iterativeWGCNA."
                        + "NOTE: assumes header for all files\n",
                        required=True)

    parser.add_argument('-c', '--classification',
                        metavar='<membership file>',
                        help="full path to gene membership file;\n"
                        + "typically 'final-membership.txt' from\n"
                        + "iterativeWGCNA output.",
                        required=True)

    parser.add_argument('-g', '--geneIdMapping',
                        metavar='<gene id mapping file>',
                        help="OPTIONAL full path to a two column, tab-delimited\n"
                        + "file providing a mapping from numeric placeholder\n"
                        + "to gene symbol or other identifier if numeric\n"
                        + "placeholders are used to identify genes in the\n"
                        + "expression dataset. Expect (ID Symbol)")

    parser.add_argument('-m', '--module',
                        metavar='<module id>',
                        help="module to extract; if not provided with\n"
                        + "return entire classified gene set\n"
                        + "use --includeUnclassified flag to return all genes")

    parser.add_argument('--addMembershipColumn',
                        help="add column to extracted expression listing\n"
                        + "module membership",
                        action='store_true')
    
    parser.add_argument('--includeUnclassified',
                        help="include unclassified genes?",
                        action='store_true')

    args = parser.parser_args()

    geneIdMap = parse_id2gene_map()
    membership = parse_membership()
    extract_expression()
    

    
