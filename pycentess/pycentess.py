from __future__ import print_function
import pycentess
import argparse
import os
import sys
from astropy.wcs import WCS
from functools import partial
import multiprocessing
import pycentess.version as pv
import pycentess.maxsect as mxs
from datetime import datetime
from tqdm import tqdm
import numpy as np






def centroid_run():
    VERSION = pv.PATHOS_VERSION
    DATEn = pv.PATHOS_UPDATED

    print('----------------------------------------------------------')
    print('pyCenTESS v.{0}'.format(VERSION))
    print('Last Update {}:'.format(DATEn))
    print('----------------------------------------------------------')


    parser = argparse.ArgumentParser(prog='pyCenTESS', description='In-/out-of transit centroid test')
    parser.add_argument('RA', type=float, nargs=1, help='Right ascension of the star')    
    parser.add_argument('DEC', type=float, nargs=1, help='Declination of the star')
    parser.add_argument('Period', type=float, nargs=1, help='Orbital period of the candidate exoplanet')
    parser.add_argument('T0', type=float, nargs=1, help='Central time of the first transit T0 (in BTJD)')
    parser.add_argument('T14', type=float, nargs=1, help='Duration of the transit (in hours)')
    parser.add_argument('min_sector', type=int, nargs=1, help='Minimum TESS Sector from which to perform the test')
    parser.add_argument('max_sector', type=int, nargs=1, help='Maximum TESS Sector from which to perform the test')


    args = parser.parse_args()

    ra = args.RA[0]
    dec= args.DEC[0]
    P  = args.Period[0]
    T0 = args.T0[0]
    T14= args.T14[0]
    min_sector = args.min_sector[0]
    max_sector = args.max_sector[0]


    if (max_sector>mxs.MAXSECT):
        print('max_sector > the maximum available sector')
        max_sector = mxs.MAXSECT
        print('max_sector= {}'.format(max_sector))

    print(' ')
    print('Starting in-/out-of transit centroid test... ')
    timestart = datetime.now()
    pycentess.stack(ra,dec,P,T0,T14,min_sector,max_sector)
