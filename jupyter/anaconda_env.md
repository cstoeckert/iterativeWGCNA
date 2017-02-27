# Sharing Anaconda Environment

## Export environment
1. Activate environment, via Anaconda Navigator or CLI
2. At the command line, execute: ```conda env export > environment.yml```
3. Share ```environment.yml``` file

## Import environment
1. Download the ```environment.yml``` file and compare it to the raw (GitHub) view to make sure your download program has not modified the file.
2. As needed, change the "name" and "prefix" in the YAML file to fit your destination directory. This is typically the first and last line, respectively. Look for your "Anaconda2/envs/" directory and change the '''/CHANGE-TO-YOUR''' part of the YAML file to fit your environment. Linux or Windows paths will work - adjust slashes as necessary.
3. At the command line, execute: ```conda env create -f environment.yml```. 
4. Activate environment, via Anaconda Navigator or CLI
5. Check the environment with: ```conda list```
