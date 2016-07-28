'''
a gene
'''

import logging

class Gene(object):
    '''
    links and manages gene properties, incuding:
    expression profile,
    module membership, 
    eigengene connectivity (kME) to assigned module
    '''

    def __init__(self, name):
        self.name = name
        self.profile = None
        self.module = None
        self.kME = None

    def one(self):
        return 1
    def two(self):
        return 2
