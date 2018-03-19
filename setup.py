from setuptools import setup, find_packages

setup(name='iterativeWGCNA',
      version='1.1.4',
      description="Iterative application of WGCNA",
      long_description='''Iterative application of
      Weighted Gene Correlation Network Analysis (WGCNA)
      to improve whole-transcriptome gene classification''',
      url='http://github.com/cstoeckert/iterativeWGCNA',
      download_url='https://github.com/cstoeckert/iterativeWGCNA/archive/v1.1.4.tar.gz',
      author='Emily Greenfest-Allen',
      author_email='allenem@pennmedicine.upenn.edu',
      license='GNU',
      packages=find_packages(),
      install_requires=['rpy2','matplotlib'],
      keywords=['network', 'WGCNA', 'gene expression', 'bioinformatics'],
      scripts=['bin/iterativeWGCNA', 'bin/iterativeWGCNA_merge'],
      zip_safe=False)
