'''
functions for defining, parsing, processing,
and validating command line arguments
'''

import argparse
from os import getcwd

def parameter_list(strValue):
    '''
    for argument parsing;
    converts a comma separated list of 'param=value' pairs
    into a parm:value hash
    '''

    params = {}
    pairs = strValue.split(',')
    for p in pairs:
        name, value = p.split('=')

        if value in ['TRUE', 'T', 'True', 't']:
            value = True
        if value in ['FALSE', 'F', 'False', 'f']:
            value = False

        params[name] = value

    return params

def restricted_float(x):
    '''
    for argument parsing; restricts float value from 0 to 1
    '''
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


def helpEpilog():
    '''
    text for help epilog
    '''

    inputFileFormatHelp = ''' 
------------------------------------------------------
Input File Format
------------------------------------------------------
iWGCNA expects a tab-delimited text file containing 
gene expression data arranged such that there is one
row per gene and one column per sample.  The first column
should contain unique gene identifiers.  For example:

GENE    Sample1    Sample2    Sample3
Gata1    500    715    1000
Phtf2    60    1000    1600
'''
    wgcnaParametersHelp = '''
------------------------------------------------------
WGCNA Parameters
------------------------------------------------------
iWGCNA can accept any parameter valid for the WGCNA
blockwiseModules function.  

See http://www.inside-r.org/packages/cran/wgcna/docs/blockwiseModules

To specify these parameters use the --wgcnaParameters flag followed by
a comma separated list of parameter=value pairs.  For example:

--wgcnaParameters maxBlockSize = 5000,corType=bicor,power=10

sets the maximum block size to 5000 genes,
the correlation type to the biweight correlation,
and the power-law scaling factor (beta) to 10

If parameters are not specified, iWGCNA uses the default WGCNA settings,
except for the following:

minModuleSize=20
saveTOMs=TRUE
minKMEtoStay=0.8
minCoreKME=0.8
reassignThreshold=0.05
networkType=signed
numericLabels=TRUE

'''

    return inputFileFormatHelp + '\n\n' + wgcnaParametersHelp


def parse_command_line_args():
    '''
    parse command line args
    '''

    parser = argparse.ArgumentParser(prog='iWGCNA',
                                     description="perform interative WGCNA analysis",
                                     epilog=helpEpilog(),
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--inputFile',
                        metavar='<gene expression file>',
                        help="full path to input gene expression file; "
                        + "if full path is not provided,\n"
                        + "assumes the file is in the working directory\n;"
                        + "see below for input file format",
                        required=True)

    parser.add_argument('-o', '--workingDir',
                        help="R working directory; where output will be saved",
                        metavar='<output dir>',
                        default=getcwd())

    parser.add_argument('-v', '--verbose',
                        help="print status messages",
                        action='store_true')

    parser.add_argument('-p', '--wgcnaParameters',
                        metavar='<param list>',
                        help="comma separated list of parameters to be passed to WGCNA's "
                        + "blockwiseModules function\n"
                        + "e.g., power=6,randomSeed=1234875\n"
                        + "see WGCNA manual & more info below",
                        type=parameter_list)

    parser.add_argument('--enableWGCNAThreads',
                        help="allow WGCNA to use threading;\nsee WGCNA manual",
                        action='store_true')

    parser.add_argument('--saveBlocks',
                        help="save WGNCA blockwise modules for each iteration",
                        action='store_true')

    parser.add_argument('--moduleMergeCutHeight',
                        help="cut height (dissimilarity threshold) for"
                        + " merging close modules after algorithm convergence;"
                        + " must be in the range (0.0, 1.0]",
                        default=0.05,
                        metavar='<dissimilarity>',
                        type=restricted_float)

    args = parser.parse_args()
    args.wgcnaParameters = set_wgcna_parameter_defaults(args.wgcnaParameters)

    return args


def set_wgcna_parameter_defaults(params):
    '''
    set default values for WGCNA blockwiseModules
    numericLabels = TRUE
    minKMEtoStay = 0.8
    '''

    if params is None:
        params = {}

    if 'numericLabels' not in params:
        params['numericLabels'] = True
    if 'minKMEtoStay' not in params:
        params['minKMEtoStay'] = 0.8
    if 'minCoreKME' not in params:
        params['minCoreKME'] = 0.8
    if 'minModuleSize' not in params:
        params['minModuleSize'] = 20
    if 'reassignThreshold' not in params:
        params['reassignThreshold'] = 0.05 # 0.0000001 # 1e-6

    return params

