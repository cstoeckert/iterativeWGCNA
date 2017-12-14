from setuptools import setup

setup(name='iterativeWGCNA',
      version='1.1.2',
      description="Iterative application of WGCNA",
      long_description='''Iterative application of
      Weighted Gene Correlation Network Analysis (WGCNA)
      to improve whole-transcriptome gene classification''',
      url='http://github.com/cstoeckert/iterativeWGCNA',
      download_url='https://github.com/cstoeckert/iterativeWGCNA/archive/v1.1.2.tar.gz',
      author='Emily Greenfest-Allen',
      author_email='allenem@pennmedicine.upenn.edu',
      license='GNU',
      packages=['iterativeWGCNA'],
      install_requires=['rpy2'],
      keywords=['network', 'WGCNA', 'gene expression', 'bioinformatics'],
      scripts=['bin/iterativeWGCNA'],
      zip_safe=False)
