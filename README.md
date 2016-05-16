# iterativeWGCNA: an WGCNA extension

## Synopsis

iterativeWGCNA provides an extension and python wrapper for the R program [Weighted Gene Co-expressionNetwork  Analysis](https://github.com/cran/WGCNA), or WGCNA, that improves the robustness of network-based classifications (modules) inferred from whole-transcriptome gene expression datasets. 

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

