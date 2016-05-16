# iterativeWGCNA: an WGCNA extension

## Synopsis

iterativeWGCNA provides an extension and Python wrapper for the R program [Weighted Gene Co-expressionNetwork  Analysis](https://github.com/cran/WGCNA), or WGCNA, that improves the robustness of network-based classifications (modules) inferred from whole-transcriptome gene expression datasets. 

## Dependencies

iterativeWGCNA has the following dependencies:

### Python (version 2.7+)

1. [rpy2](https://pypi.python.org/pypi/rpy2): a Python interface for R (v. 2.7.9+)
  * will be installed by Python setup tools or package manager if the iwgcna packages is installed

### R language for statistical computing

[R](https://cran.r-project.org/) version 3.* must be available on the system and the binary executable in the system PATH

1. [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall): Weighted Gene Co-expression Network Analysis package and Bioconductor dependencies

## Code Example

iterativeWGCNA can be run without installing the iwgcna package by excuting the wrapper script `run_iwgcna.py` in the iterativeWGCNA directory. At a minimum, the `-i` option (`--inputFile`) denoting the full path to the input file must be specified.

```
python run_iwgcna.py -i <input_file_path> 
```

Execute `run_iwgcna.py` with the `-h` (`--help`) flag to see all command line options and file format descriptions:

```
python run_iwgcna.py -h
```

if the iwgcna package was installed, iterativeWGCNA can also be run at the package level using the `-m` switch:

```
python -m iwgcna -h
```

## Installation

iterativeWGCNA can be run without installing the iwgcna package.  To install the package, use the following command:

```
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install
```

