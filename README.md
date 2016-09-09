# iterativeWGCNA: a WGCNA extension

## Synopsis

iterativeWGCNA provides an extension and Python wrapper for the R program [Weighted Gene Co-expression Network Analysis](https://github.com/cran/WGCNA), or WGCNA, that improves the robustness of network-based classifications (modules) inferred from whole-transcriptome gene expression datasets.

## Dependencies

iterativeWGCNA has the following dependencies:

### Python (version 2.7+)

1. [rpy2](https://pypi.python.org/pypi/rpy2): a Python interface for R (v. 2.7.9+)
  * will be installed by Python setup tools or package manager if the iterativeWGCNA packages is installed

### R language for statistical computing

[R](https://cran.r-project.org/) version 3.* must be available on the system and the binary executable in the system PATH. Note that if you are installed WGCNA from CRAN, the most recent version of R you can run (that supports WGCNA) is 3.2.0.

1. [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall): Weighted Gene Co-expression Network Analysis package and Bioconductor dependencies

## Code Example

iterativeWGCNA can be run without installing the iterativeWGCNA package by excuting the wrapper script `run_iterative_wgcna.py` in the iterativeWGCNA directory. At a minimum, the `-i` option (`--inputFile`) denoting the full path to the input file must be specified.

```
python run_iterative_wgcna.py -i <input_file_path> 
```

Execute `run_iterative_wgcna.py` with the `-h` (`--help`) flag to see all command line options and file format descriptions:

```
python run_iterative_wgcna.py -h
```

if the iterativeWGCNA package was installed, iterativeWGCNA can also be run at the package level using the `-m` switch:

```
python -m iterativeWGCNA -h
```

## Installation

### Instructions

iterativeWGCNA can be run without installing the iterativeWGCNA package.  To install the package, use the following command:

```
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install
```

To install iterativeWGCNA in a user directory (non-system), then the `--user` switch will install the package to `.local/bin` (on Linux systems):

```
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install --user
```

### Workarounds

Access to the `readline` library in the context of the `rpy2` library can be problematic and has been [reported elsewhere](https://github.com/ContinuumIO/anaconda-issues/issues/152). In trying to run iterativeWGCNA, an error like the following would be observed:
```
Traceback (most recent call last):
  File "/home/_USER_/iterativeWGCNA-master/run_iterative_wgcna.py", line 7, in <module>
    from iterativeWGCNA.iterativeWGCNA import IterativeWGCNA
  File "/home/_USER_/iterativeWGCNA-master/iterativeWGCNA/iterativeWGCNA.py", line 17, in <module>
    import rpy2.robjects as ro
  File "/home/_USER_/bin/anaconda2/lib/python2.7/site-packages/rpy2/robjects/__init__.py", line 15, in <module>
    import rpy2.rinterface as rinterface
  File "/home/_USER_/bin/anaconda2/lib/python2.7/site-packages/rpy2/rinterface/__init__.py", line 100, in <module>
    from rpy2.rinterface._rinterface import *
ImportError: /home/_USER_/bin/anaconda2/lib/python2.7/site-packages/rpy2/rinterface/../../../../libreadline.so.6: undefined symbol: PC
```

The workaround is to uncomment the following line in `run_iterative_wgcna.py`:
```
# import readline
```
