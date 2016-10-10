# pylint: disable=invalid-name
'''
manage colors
'''

import logging
from random import randint
from matplotlib import colors as ref_colors


class Colors(object):
    '''
    manage colors for modules
    '''
    def __init__(self):
        # base colors taken from WGCNA standard coloring for modules
        self.logger = logging.getLogger('iterativeWGCNA.Colors')

        self.base_colors = ['turquoise', 'blue', 'brown', 'yellow', 'green',
                            'red', 'black', 'pink', 'magenta', 'purple',
                            'greenyellow', 'tan', 'salmon', 'cyan', 'midnightblue',
                            'lightcyan', 'lightgreen', 'lightyellow', 'royalblue',
                            'darkred', 'darkgreen', 'darkturquoise', 'orange',
                            'darkorange', 'skyblue', 'saddlebrown', 'steelblue',
                            'paleturquoise', 'violet', 'darkolivegreen',
                            'darkmagenta']

        self.used_colors = []


    def assign_color(self, n):
        '''
        assigns a color
        if n <= len(base_colors) assigns the base color whose index is n - 1
        else generates a random color
        '''
        color = None
        if n <= len(self.base_colors):
            color = ref_colors.cnames[self.base_colors[n - 1]] # get hex representation
        else:
            color = self.__generate_random_color()

        self.used_colors.append(color)
    
        return color


    def __generate_random_color(self):
        '''
        generate a random color
        '''
        color = '#' + '%06X' % randint(0, 0xFFFFFF)
        while color in self.used_colors:
            color = '#' + '%06X' % randint(0, 0xFFFFFF)
        return color
