#!/usr/bin/env python

"""
Jean coupon - 2018
script to load dust map and get extinction value
"""

from __future__ import print_function, division

import os, sys
import numpy as np

from astropy.io import ascii,fits
import astropy.wcs as wcs

from astropy import units as u
from astropy.coordinates import SkyCoord


#np.random.seed(20091982)


"""

-------------------------------------------------------------
main functions
-------------------------------------------------------------

"""

def main(args):
    """ Routine to get and apply dust correction
    """

    # column names
    ra_name, dec_name = args.columns.split(',')

    # load dust maps
    sMap, nMap, sWcs, nWcs = getMaps(args)

    # input catalogue
    fileInName = args.input
    fileIn = fits.open(fileInName)
    data = fileIn[1].data

    N = len(data)

    # ra = data['RA']
    # dec = data['DEC']


    EB_V = np.zeros(N)
    x = np.zeros(N)
    y = np.zeros(N)

    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')

    north = coord.galactic.b.degree > 0.0
    south = coord.galactic.b.degree < 0.0

    if len(coord[north]) > 0:
        x[north], y[north]  = wcs.utils.skycoord_to_pixel(
            coord[north], nWcs,  origin=0)
    if len(coord[south]) > 0:
        x[south], y[south]  = wcs.utils.skycoord_to_pixel(
            coord[south], sWcs,  origin=0)

    for i in range(N):
        if north[i]:
            EB_V[i] = nMap[iround(y[i]), iround(x[i])]
        if south[i]:
            EB_V[i] = sMap[iround(y[i]), iround(x[i])]

    if args.correct:
        sys.stderr.write("Correcting fluxes...\n")
        correctBands(args, data, EB_V)

    cols = []
    cols.append(fits.Column(name='EB_V', format='E', array=EB_V))

    if args.add_corr:
        sys.stderr.write("Add correction columns...\n")
        addCorr(args, data, EB_V, cols)

    tbhdu = fits.BinTableHDU.from_columns(
        fileIn[1].columns + fits.ColDefs(cols))
    tbhdu.writeto(args.output, clobber=True)

    return



def correctBands(args, data, EB_V):

    band = str.split(args.band,    ",")
    coef = str.split(args.coef,    ",")
    N = len(band)

    if N == 0: return

    if len(band) != len(coef):
        raise ValueError("band and coef must have the same number of elements")

    for i in range(N):
        data[band[i]] *=  pow(10.0, +0.4*float(coef[i]) * EB_V)


    return

def addCorr(args, data, EB_V, cols):

    band = str.split(args.band,    ",")
    coef = str.split(args.coef,    ",")
    N = len(band)

    if N == 0: return

    if len(band) != len(coef):
        raise ValueError("band and coef must have the same number of elements")

    for i in range(N):
        corr =  pow(10.0, +0.4*float(coef[i]) * EB_V)
        cols.append(
            fits.Column(name='EB_V_corr_'+band[i], format='E', array=corr))

    return


"""

-------------------------------------------------------------
utils
-------------------------------------------------------------

"""

def getMaps(args):

    sFile = fits.open(args.southFile)
    nFile = fits.open(args.northFile)

    sMap = sFile[0].data
    nMap = nFile[0].data

    sWcs = wcs.WCS(sFile[0].header)
    nWcs = wcs.WCS(nFile[0].header)

    return sMap, nMap, sWcs, nWcs

def iround(x):
    """
    iround(number) -> integer
    Round a number to the nearest integer.
    From https://www.daniweb.com/software-development/python/threads/299459/round-to-nearest-integer-
    """

    return int(round(x) - .5) + (x > 0)

"""

-------------------------------------------------------------
main
-------------------------------------------------------------

"""



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help="input file")
    parser.add_argument('output', help="output file")

    parser.add_argument(
        '-c', '--columns',
        default="ra,dec",
        help='Names of coordinates colums (default: ra,dec)')


    parser.add_argument(
        '-s', '--southFile',
        default="/Users/coupon/data/SchlegelDust/SFD_dust_4096_sgp.fits",
        help='Dust map (south galactic)')
    parser.add_argument(
        '-n', '--northFile',
        default="/Users/coupon/data/SchlegelDust/SFD_dust_4096_ngp.fits",
        help='Dust map (north galactic)')

    parser.add_argument(
        '-band',
        help="List of bands to correct the extinction for")
    parser.add_argument(
        '-coef',
        help="Corresponding Albda/E(B-V) (in mags)")

    parser.add_argument(
        '-correct', action='store_true',
        help='Correct fluxes for extinction')
    parser.add_argument(
        '-add_corr', action='store_true',
        help='Add flux correction factor for extinction')

    args = parser.parse_args()

    main(args)
