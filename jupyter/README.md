# Interactive Demonstration of iterativeWGCNA and WGCNA
____

Available are [Jupyter](http://jupyter.org/) notebooks that will allow you to run iterativeWGCNA and WGCNA interactively in your web browser on your __local__ machine.

To reduce your pain in setting up a suitable environment, we recommend that you use [Anaconda2](https://www.continuum.io/downloads) to run Jupyter and notebooks. We have provided *environment files* that should recapitulate our system onto yours, for Linux and Windows. You could also try Miniconda or go the manual route.

## iterativeWGNCA
Details to come... [try notebook](iterativeWGCNA.ipynb)

## WGNCA
Details to come... [try notebook](Standard_WGCNA.ipynb)

## Anaconda

### Importing the environment
1. Download the ```anaconda2_env_***.yml``` file suitable to your operating system ([windows](anaconda2_env_win.yml) vs [linux](anaconda2_env_lin.yml)) and compare it to the raw (GitHub) view to make sure your download program has not modified the file.
2. As needed, change the "prefix" in the YAML file to fit your destination directory. This is typically the last line. Look for your "Anaconda2/envs/" directory and change the '''/CHANGE-TO-YOUR''' part of the YAML file to fit your environment and system.
3. At the command line, execute: ```conda env create -f anaconda2_env.yml```.
4. Activate environment, via Anaconda Navigator or CLI
5. Check the environment with: ```conda list```

### Exporting the environment
(This is mainly for us to remember how to do this!)
1. Activate environment, via Anaconda Navigator or CLI
2. At the command line, execute: ```conda env export -n iterativeWGCNA > anaconda2_env.yml```
3. Share ```environment.yml``` file
