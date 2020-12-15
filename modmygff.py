#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
from gff3 import Gff3

from pprint import pprint

eg_gff_path = os.path.join(os.getcwd(), 'data', 'Slin_CCMP2456',
                           'S.linucheae_CCMP2456_eg1.gff')

eg_gff: Gff3 = Gff3(gff_file=eg_gff_path)
# pprint(eg_gff.lines[0])

eg_gff_ext_path = os.path.join(os.getcwd(), 'data', 'Slin_CCMP2456',
                               'S.linucheae_CCMP2456_eg1_ext.gff')

with open(eg_gff_ext_path, 'w') as eg_gff_ext_file:
    eg_gff.write(eg_gff_ext_file)
