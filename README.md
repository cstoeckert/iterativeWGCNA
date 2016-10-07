# iterativeWGCNA: a WGCNA extension

## Synopsis

iterativeWGCNA provides an extension and Python wrapper for the R program [Weighted Gene Co-expression Network Analysis](https://github.com/cran/WGCNA), or WGCNA, that improves the robustness of network-based classifications (modules) inferred from whole-transcriptome gene expression datasets.

## Contents

1. [Setup and Installation](#setup)
2. [Usage](#usage)
A. [Run iterativeWGCNA](#run)
B. [Summarize results](#summarize)
3. [Troubleshooting](#troubleshooting)

<a name="setup"></a>
## Setup and Installation

### A. Dependencies

iterativeWGCNA has the following dependencies:

#### i. Python (version 2.7+)* 

1. [rpy2](https://pypi.python.org/pypi/rpy2): a Python interface for R (v. 2.7.9+)
2. [matplotlib](http://matplotlib.org/): a python 2D plotting library

* will be installed by Python setup tools or package manager if the iterativeWGCNA packages is installed

#### ii. R language for statistical computing

[R](https://cran.r-project.org/) version 3.* must be available on the system and the binary executable in the system PATH. Note that if you installed WGCNA from CRAN, the most recent version of R you can run (that supports WGCNA) is 3.2.0.

1. [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall): Weighted Gene Co-expression Network Analysis package and Bioconductor dependencies
2. [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html): Pretty Heatmaps

### B. Installation

iterativeWGCNA can be run without installing the package as long as the requisite Python dependencies (rpy2 and matplotlib) are already present on the system.  Installing the package will also install the dependencies.  To install the iterativeWGCNA package run the following command:

```
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install
```

NOTE: depending on your system this may require administrative (e.g., sudo) permissions.  As a work around, specify the `--user` switch to install iterativeWGCNA and its dependencies to a local (user) library (e.g., `.local/bin` on a Linux system) as follows:

```
git clone https://github.com/cstoeckert/iterativeWGCNA.git
cd iterativeWGCNA
python setup.py install --user
```


## 

### [How to run iterativeWGCNA](docs/basic_usage.md)
### [Summarize iterativeWGCNA output](docs/summaries.md)


## Installation


<a name="troubleshooting"></a>
### Troubleshooting

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
