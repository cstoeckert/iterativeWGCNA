### 1.1.6 / 2018-03-19
  * ability to gzip TOM .RData files added `--gzipTOMs` option
  * new release on PyPI 
  * bug fix: saveTOMs disabled by default; documentation updated
  * bug fix: issue parsing boolean WGCNA parameters (e.g. saveTOMs=FALSE or cosineCorrelation=TRUE) resolved 
  
###  1.1.3
  * added script to adjust final module merge
	* see [Add-ons](/README.md#add-ons) and updated [Output Files](/README.md#output-files) for more information
  * fixed Python 3.3+ bug with converting odict_values to ro.StrVector
  * added `--debug` option; currently only prints extensive debugging statements for module merge stage

### 1.1.0 / 2017-06-28
* Fix final module merge to improve efficiency and recalculate eigengenes after each merge 
* Change program output to generate iteration-specific output files in easily accessible directory structure

### 1.0.0 / 2017-02-06
* First numbered release
* Version names will follow "Semantic Versioning" of X.Y.Z where:
    * Increment Z when you fix something
    * Increment Y when you add a new feature
    * Increment X when you break backwards-compatibility or add major features
* Tags/Releases will be against Master, as a general guideline to keep it simple
