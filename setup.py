from setuptools import setup

setup(name='iterativeWGCNA',
      version='0.1',
      description="Iterative application of WGCNA",
      long_description='''Iterative application of
      Weighted Gene Coexpression Network Analysis (WGCNA)
      to improve whole-transcriptome gene classification''',
      url='http://github.com/cstoeckert/iterativeWGCNA',
      author='Emily Greenfest-Allen',
      author_email='allenem@upenn.edu',
      license='GNU',
      packages=['iterativeWGCNA'],
      install_requires=['rpy2'],
      zip_safe=False)
