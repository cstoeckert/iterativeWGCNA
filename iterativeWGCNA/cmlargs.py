#pylint: disable=anomalous-backslash-in-string
#pylint: disable=invalid-name

'''
functions for defining, parsing, processing,
and validating command line arguments
'''

import re
import argparse
from os import getcwd
from .io.utils import warning

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

        # Test and cast for booleans
        if value.upper() in ['TRUE', 'T']:
            params[name] = True
        elif value.upper() in ['FALSE', 'F']:
            params[name] = False
        # Test and cast for integer
        elif value.isdigit():
            params[name] = int(value)
        # Test and cast for float
        elif re.match("^\d+?\.\d+?$", value):
            params[name] = float(value)
        else:
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


def summaryHelpEpilog():
    '''
    text for help epilog for
    summary
    '''
    # TODO Help for summarize script
    return "COMING SOON"

def helpEpilog():
    '''
    text for help epilog
    '''

    inputFileFormatHelp = '''
------------------------------------------------------
Input File Format
------------------------------------------------------
iterativeWGCNA expects a tab-delimited text file containing
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
iterativeWGCNA can accept any parameter valid for the WGCNA
blockwiseModules function.

See http://www.inside-r.org/packages/cran/wgcna/docs/blockwiseModules

To specify these parameters use the --wgcnaParameters flag followed by
a comma separated list of parameter=value pairs.  For example:

--wgcnaParameters maxBlockSize=5000,corType=bicor,power=10

sets the maximum block size to 5000 genes,
the correlation type to the biweight correlation,
and the power-law scaling factor (beta) to 10

If parameters are not specified, iterativeWGCNA uses the default WGCNA settings,
except for the following:

minModuleSize=20
saveTOMs=TRUE
minKMEtoStay=0.8
minCoreKME=0.8
networkType=signed
numericLabels=TRUE

'''

    return inputFileFormatHelp + '\n\n' + wgcnaParametersHelp


def parse_command_line_args(program='iterativeWGCNA', description='perform iterativeWGCNA analysis'):
    '''
    parse command line args
    '''

    parser = argparse.ArgumentParser(prog=program,
                                     description=description,
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

    parser.add_argument('--debug',
                        help="print debugging messages",
                        action='store_true')
    
    parser.add_argument('-p', '--wgcnaParameters',
                        metavar='<param list>',
                        help="comma separated list of parameters to be passed to WGCNA's "
                        + "blockwiseModules function\n"
                        + "e.g., power=6,randomSeed=1234875\n"
                        + "see WGCNA manual & more info below",
                        type=parameter_list)

    parser.add_argument('--enableWGCNAThreads',
                        help="enable WGCNA to use threading;\nsee WGCNA manual",
                        action='store_true')

    parser.add_argument('--skipSaveBlocks',
                        help="do not save WGCNA blockwise modules for each iteration;\n"
                        + "NOTE: without blocks summary graphics cannot be generated.\n"
                        + "Also will not saveTOMs.",
                        action='store_true')

    parser.add_argument('-f', '--finalMergeCutHeight',
                        help="cut height for final merge (after iterations are assembled)",
                        default=0.05,
                        metavar='<cut height>',
                        type=restricted_float)

    args = parser.parse_args()
    args.wgcnaParameters = set_wgcna_parameter_defaults(args.wgcnaParameters, args.skipSaveBlocks)

    return args


def set_wgcna_parameter_defaults(params, skipSaveBlocks):
    '''
    set default values for WGCNA blockwiseModules
    numericLabels = TRUE
    minKMEtoStay = 0.8
    '''

    if params is None:
        params = {}

    params['numericLabels'] = True # override user choice
    
    if 'networkType' not in params:
        params['networkType'] = 'signed'
    if 'minKMEtoStay' not in params:
        params['minKMEtoStay'] = 0.8
    if 'minCoreKME' not in params:
        params['minCoreKME'] = params['minKMEtoStay']
    if 'minModuleSize' not in params:
        params['minModuleSize'] = 20
    if 'reassignThreshold' not in params:
        params['reassignThreshold'] = 0.05 # 0.000001 # 1e-6
    if 'power' not in params:
        params['power'] = 6
    if 'saveTOMs' not in params:
        params['saveTOMs'] = True
    if skipSaveBlocks: # if blocks are not saved; TOMs are not saved
        params['saveTOMs'] = False

    return params


def parse_summary_command_line_args():
    '''
    parse command line args for summary
    '''

    parser = argparse.ArgumentParser(prog='iterativeWGCNA network summary',
                                     description="generate graphical results from interative WGCNA analysis",
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

    parser.add_argument('-p', '--power',
                        metavar='<power law beta>',
                        help="power law beta for weighting the adjacency matrix",
                        default=6,
                        type=int)

    parser.add_argument('--signed',
                        help="generate signed adjacency matrix?",
                        action='store_true')

    parser.add_argument('--minKMEtoStay',
                        help="provide minKMEtoStay used for network generation",
                        default=0.80,
                        metavar='<minKMEtoStay>',
                        type=restricted_float)

    parser.add_argument('--enableWGCNAThreads',
                        help="enable WGCNA to use threading;\nsee WGCNA manual",
                        action='store_true')

    parser.add_argument('--generateNetworkSummary',
                        metavar='<view type>',
                        choices=['all', 'network', 'input'],
                        help="generate summary overview of the network (dendrogram & heatmap):\n"
                        + "network - network comprised only of classified genes\n"
                        + "input - all genes, with classified highlighted by module assignments\n"
                        + "all - both views\n"
                        + "NOTE: all adjacency matrix calculations are\n"
                        + "done in one block and may fail due to memory allocation\n"
                        + "issues for large gene-sets")

    parser.add_argument('-e', '--edgeWeight',
                        metavar='<min edge weight>',
                        default=0.5,
                        help="min edge weight for network summary; filters for\n"
                        + "connections supported by a correlation >= threshold",
                        type=restricted_float)


    return parser.parse_args()


