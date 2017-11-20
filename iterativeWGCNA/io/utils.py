'''
I/O Utils
'''
from __future__ import print_function
from __future__ import with_statement

from sys import stderr
from math import isnan
import os
import re
import rpy2.robjects as ro
from ..r.imports import rsnippets


def xstr(value):
    '''
    handle nulls/nan in string conversion
    '''
    if value is None:
        return ''
    if value == 'NULL':
        return ''
    if isnan(value):
        return 'NA'

    return str(value)


def warning(*objs):
    '''
    wrapper for writing to stderr
    '''
    print(*objs, file=stderr)
    stderr.flush()


def create_dir(dirName):
    '''
    check if directory exists in the path, if not create
    '''
    try:
        os.stat(dirName)
    except OSError:
        os.mkdir(dirName)

    return dirName


def write_data_frame(df, fileName, rowLabel):
    '''
    write data frame to file; creates new file
    if none exists, otherwise appends new eigengenes
    to existing file
    '''
    try:
        os.stat(fileName)
    except OSError:
        header = (rowLabel,) + tuple(df.colnames)
        with open(fileName, 'w') as f:
            print('\t'.join(header), file=f)
    finally:
        df.to_csvfile(fileName, quote=False, sep='\t', col_names=False, append=True)

        
def read_data(fileName):
    '''
    read gene expression data into a data frame
    and convert numeric (integer) data to real
    '''
    data = ro.DataFrame.from_csvfile(fileName, sep='\t', header=True, row_names=1)
    return rsnippets.numeric2real(data)


def transpose_file_contents(fileName, rowLabel):
    '''
    read in a file to a dataframe, transpose, and output
    use this instead of R transforms b/c R will concatenate
    an "X" to gene ids starting with a number
    '''
    with open(fileName, 'r') as f:
        content = [line.rstrip().split() for line in f]

    header = True
    with open(fileName, 'w') as f:
        for line in zip(*content):
            if header:
                line = list(line)
                line[0] = rowLabel
                line = tuple(line)
                header = False

            # b/c rpy2 replaces '-' in gene symbols with '.'
            if '.' in line[0]:
                line = list(line)
                line[0] = line[0].replace('.', '-')
                line = tuple(line)

            # b/c R tacks an X on to gene names that start
            # with a #
            if re.search('^X\d', line[0]) is not None:
                print(line[0], file=stderr)
                line = list(line)
                line[0] = re.sub('^X', '', line[0])
                line = tuple(line)

            print('\t'.join(line), file=f)
