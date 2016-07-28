# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
# pylint: disable=unused-import

'''
Gene properties + classification check
'''

import logging # for debugging purposes

class (object):
    '''
    links and manages gene properties,
    incuding:
      expression profile,
      module membership,
      eigengene connectivity (kME) to assigned module
    '''

    def __init__(self, profile):
        self.profile = profile
        self.module = None
        self.kME = None

    def is_classified(self):
        '''
        check if gene is classified
        '''
        return self.module != 'UNCLASSIFIED'
