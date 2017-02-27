# Sharing Anaconda Environment

## Export environment
1. Activate environment, via Anaconda Navigator or CLI
2. At the command line, execute: ```conda env export > environment.yml```
3. Share ```environment.yml``` file

## Import environment
1. At the command line, execute: ```conda env create -f environment.yml```
2. Activate environment, via Anaconda Navigator or CLI
3. Check the environment with: ```conda list```
