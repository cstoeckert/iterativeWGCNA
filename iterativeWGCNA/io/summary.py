'''
progress summary output
'''
from __future__ import print_function
import os

def write_gene_counts(iteration, initial, fit):
    '''
    writes the number of kept and dropped genes at the end of an iteration
    '''
    fileName = 'gene-counts.txt'
    try:
        os.stat(fileName)
    except OSError:
        header = ('Iteration', 'Initial', 'Fit', 'Residual')
        with open(fileName, 'a') as f:
            print('\t'.join(header), file=f)
    finally:
        with open(fileName, 'a') as f:
            print('\t'.join((iteration, str(initial),
                             str(fit), str(initial - fit))), file=f)

