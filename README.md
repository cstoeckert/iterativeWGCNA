# iterativeWGCNA: a WGCNA extension

## Synopsis

iterativeWGCNA provides a Python-wrapped extension for the R program [Weighted Gene Correlation Network Analysis](https://github.com/cran/WGCNA) (WGCNA) that improves the robustness of network-based classifications (modules) inferred from whole-transcriptome gene expression datasets.

## Contents


### Setup and Installation

* [Dependencies](#dependencies)
* [Installation](#installation)
  
### Usage

* [Running iterativeWGCNA](#running-iterativewgcna)
  
### Troubleshooting

* [libreadline.so.6: undefined symbol](#libreadlineso6-undefined-symbol)
* [Cannot install rpy2 on OSX](#cannot-install-rpy2-with-latest-r-version-34x-on-macos)

## Setup and Installation

### Dependencies

iterativeWGCNA has the following dependencies:

#### R language for statistical computing

[R](https://cran.r-project.org/) version 3.* must be available on the system and the binary executable in the system PATH.

> NOTE: the most recent version of R that supports WGCNA is 3.3.x

iterativeWGCNA requires that the following R packages be installed:

* [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall): Weighted Gene Co-expression Network Analysis package and Bioconductor dependencies

#### Python

iterativeWGCNA requires Python version 2.7 or higher.  It is designed to be future compatible with Python 3+.  iterativeWGCNA requires the following Python packages:

* [rpy2](https://pypi.python.org/pypi/rpy2): a Python interface for R (v. 2.7.9+)

> NOTE: the most recent version of rpy2 requires python 3.x

If missing, rpy2 will be installed by the iterativeWGCNA installer.  See below.

### Installation

iterativeWGCNA can be run without installing the package as long as the requisite Python dependencies (rpy2) and R are already present on the system.  Installing the package will install any missing *Python* dependencies.

> iterativeWGCNA is reposited in the Python Package Index (PyPI) and can be installed via `pip` or `easy_install`.

To install the iterativeWGCNA package from an archived release or from git, use the `python setup.py` script as folows:

```bash
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install
```

> NOTE: depending on your system this may require administrative (e.g., sudo) permissions.

As a work around, specify the `--user` switch to install iterativeWGCNA and its dependencies to a local (user) library (e.g., `.local/bin` on a Linux system) as follows:

```sh
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install --user
```

## Usage

### Running iterativeWGCNA

1. [Quick Start](#quick-start)
1. [Command Line Options](#command-line-options)
1. [WGCNA Parameters](#wgcna-parameters)
1. [Input File Format](#input-file-format)
1. [Output Files](#output-files)

#### Quick Start

iterativeWGCNA can be run without installing the iterativeWGCNA package by excuting the wrapper script `run_iterative_wgcna.py` in the iterativeWGCNA directory. At a minimum, the `-i` option (`--inputFile`) denoting the full path to the input file must be specified.

```sh
python run_iterative_wgcna.py -i <input_file_path> 
```

if the iterativeWGCNA package was installed, iterativeWGCNA can also be run at the package level using the `-m` switch:

```sh
python -m iterativeWGCNA -i <input_file_path> 
```

#### Command Line Options

Execute `run_iterative_wgcna.py` with the `-h` (`--help`) switch to see all command line options and additional usage information, including details on file formats.

```sh
python run_iterative_wgcna.py -h
```

```diff
-h, --help
   show help message and exit
    
-i <gene expression file>, --inputFile <gene expression file>
   full path to input gene expression file; if full path is not provided,
   assumes the file is in the working (output) directory
+ required
   
-o <output dir>, --workingDir <output dir>
   R working directory; where output will be saved
   default: current directory
   
-v, --verbose
   print status messages

-p <param list>, --wgcnaParameters <param list>
   comma separated list of parameters to be passed to WGCNA's blockwiseModules function
   e.g., power=6,randomSeed=1234875
   see 'blockwiseModules' section of the WGCNA manual for more information
   
--enableWGCNAThreads
    enable WGCNA to use threads
    
--skipSaveBlocks
    do not save WGCNA blockwise modules for each iteration
	NOTE: blocks are necessary to generate summary graphics
	
-f, --finalMergeCutHeight <cut height>
	cut height (max dissimilarity) for final module merge
	(after algorithm convergence); [0, 1.0], default=0.05
    
```

#### WGCNA Parameters

iterativeWGCNA can accept any parameter valid for the WGCNA blockwiseModules function.  See https://www.rdocumentation.org/packages/WGCNA/versions/1.41-1/topics/blockwiseModules for full details

> To specify these parameters use the `--wgcnaParameters` flag followed by a comma separated list of parameter=value pairs.

For example:

`--wgcnaParameters maxBlockSize=5000,corType=bicor,power=10`

sets the maximum block size to 5000 genes, the correlation type to the biweight correlation, and the power-law scaling factor (beta) to 10

> WGCNA's `blockwiseModules` function partitions the gene set into a set of blocks each containing at most `maxBlockSize` genes.

*To run iterativeWGCNA in a single block, set `maxBlockSize` to a value > than the number of genes in your geneset*.

> NOTE: for large datasets (>10,000 genes or probes), adjacency and TOM matrix calculations done in a single block may fail due to memory allocation issues 

see the [WGCNA large dataset tutorial, section 2.c.2](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf) for more details

If WGCNA parameters are not specified, iterativeWGCNA uses the default WGCNA settings for the `blockwiseModules` function, except for the following:

```python
minModuleSize = 20 # minimum number of genes in a detected module
saveTOMs = TRUE # save the topological overlap matrices for each block in the block data structure
minKMEtoStay = 0.8 # minimum eigengene connectivity (kME) required for a gene to be retained in its assigned module
minCoreKME = 0.8 # if the module does not have minModuleSize genes with eigengene connectivity at least minCoreKME, the module is disbanded
reassignThreshold = 0.05 # per-set p-value ratio threshold for reassigning genes between modules
networkType = 'signed' # character string specifying network type. Allowed values are "unsigned", "signed", and "signed hybrid"
numericLabels = TRUE # label modules by numbers (e.g., 0,1,2) instead of colors
```

#### Input File Format

iterativeWGCNA expects a `tab-delimited` text file containing gene expression data arranged such that there is one row per gene and one column per sample.  The first column should contain `unique` gene identifiers.  For example:

| Gene | Sample1 | Sample2 | Sample3 |
| --- | --- | --- | --- |
| Gata1 | 500 | 715 | 1000 |
| Phtf2 | 60 | 1000 | 1600 |

> NOTE: We recommend using numeric gene identifiers to uniquely label genes in the input file as R will do some character substitutions (e.g., '.' for '-') and concatenate 'X' to gene symbols starting with a number, leading to erroneous mapping between data frames and potential loss of data.


> iterativeWGCNA will accept `gzipped` input files.


#### Output Files

An **iteration** of iterativeWGCNA comprises one run of blockwiseWGCNA followed by an eigengene-connectivity (kME) goodness of fit assessment.  A **pass** of iterativeWGCNA comprises multiple iterations applied to an expression dataset until no more residuals to the kME-fit are found.  A new pass is initiated by creating a new expression dataset from all residuals to the kME-fit found during the previous pass.

> Modules are uniquely identified by the numerical assignment and the iteration in which they were first detected: e.g., `P1_I2_M1` is module 1, detected in the second iteration of the first pass.  **Unclassified genes are labeled UNCLASSIFIED with a kME of NA**.

Results from each pass and iteration are saved in a series of directories, labeled as:

> passN: results from the numbered (N) pass
> iN: results from the numbered (N) iteration

The directory structure and output files are as follows:

```
├── output_directory
│   ├── iterativeWGCNA.log: main log file for the iterativeWGCNA run
│   ├── iterativeWGCNA-R.log: log file for R; catches R errors and R warning messages
│   ├── gene-counts.txt: tally of number of genes fit and residual to the fit with each iteration
│   ├── final-eigengenes.txt: eigengenes for final modules after final network assembly and merging of close modules
│   ├── final-kme-histogram.pdf: histogram of eigengene connectivities (kME) in the final classification
│   ├── final-membership.txt: gene-module assignments and kME after final iteration and merge
│   ├── passN
│   │   ├── initial-pass-expression-set.txt: pass input
│   │   ├── kme_histogram.pdf: histogram of eigengene connectivities for genes classified during pass
│   │   ├── membership.txt: gene-module assignments and kME for genes classfied during pass
│   │   ├── iN
│   │   │   ├── eigengenes.txt: eigengenes of modules detected during the iteration
│   │   │   ├── kme_histogram.pdf: kME histogram after pruning of WGCNA result based on kME
│   │   │   ├── membership.txt: gene membership after kME-based goodness of fit (Pruning)
|   │   │   ├── summary.txt: summaries pass (number genes input, classfied, residual, and number of detected modules)
│   │   |   ├── wgcna-blocks.RData: R data object containing input expression data (expression) and results from blockwise WGCNA (blocks)
│   │   │   ├── wgcna-kme_histogram.pdf: kME histogram based on WGCNA classification
│   │   │   ├── wgcna-membership.txt: gene membership from WGCNA classification
```

## Troubleshooting

### libreadline.so.6 undefined symbol

Access to the `readline` library in the context of the `rpy2` library and an Ananconda install can be problematic and has been [reported elsewhere](https://github.com/ContinuumIO/anaconda-issues/issues/152). In trying to run iterativeWGCNA, an error like the following would be observed:
```
Traceback (most recent call last):
  File "../iterativeWGCNA-master/run_iterative_wgcna.py", line 7, in <module>
    from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA
  File "../iterativeWGCNA-master/iterativeWGCNA/iterativeWGCNA.py", line 17, in <module>
    import rpy2.robjects as ro
  File "../lib/python2.7/site-packages/rpy2/robjects/__init__.py", line 15, in <module>
    import rpy2.rinterface as rinterface
  File "../lib/python2.7/site-packages/rpy2/rinterface/__init__.py", line 100, in <module>
    from rpy2.rinterface._rinterface import *
ImportError: ../lib/python2.7/site-packages/rpy2/rinterface/../../../../libreadline.so.6: undefined symbol: PC
```

The workaround is to uncomment the readline import in the `run_iterative_wgcna.py` script:

```python
# import readline
```

### Cannot install rpy2 with latest R (version 3.4.x) on macOS

#### error: command 'gcc' failed with exit status 1

The build process for many R libraries explicitly uses ```gcc```, instead of the system default, which on OSX is ```clang```.  To override the default compiler, set the CC environmental variable as folows:

```bash
export CC=clang
```

#### clang: error: unsupported option '-fopenmp'

This is a known issue with an open ticket in the rpy2 project (see issue [#403](https://bitbucket.org/rpy2/rpy2/issues/403/cannot-pip-install-rpy2-with-latest-r-340)).  R 3.4.0 was built using the ```-fopenmp``` flag with Clang 4.0.0, which is not supplied by Apple.  There are several suggested workarounds (e.g., installing the LLVM library via homebrew) that do not work for all system configurations.  We recommend [downloading the rpy2 source](https://bitbucket.org/rpy2/rpy2/src), unpacking, and editing the ```setup.py``` file after line 268 (the line above the comment ```# OS X's frameworks need special attention```) as follows:

```python
  if "-fopenmp" in unknown:  # remove linker argument
        unknown.remove("-fopenmp")
```

With this fix you should be able to build rpy2 from the downloaded source as follows:

```bash 
python setup.py install
```




