# import libraries
import logging as logger
import numpy as np
from collections import defaultdict
import logging
import os
import csv


import sys
sys.path.insert(0, '/master/')
sys.path.insert(0, '/master/lib')
sys.path.insert(0, '/master/lines')

# TODO: FIX THESE IMPORT ISSUES
import master.main
import master.lib
import master.lines


lgr = logging.getLogger('Default Debug Logger')
lgr.setLevel(logging.DEBUG)

fh = logging.FileHandler('debug_log.log')
fh.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)

# add the handlers to the logger
lgr.addHandler(fh)
lgr.addHandler(ch)

