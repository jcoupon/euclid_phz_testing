#!/usr/bin/env python

"""
Jean coupon - 2018
scripts to simulate photometry
"""

import os, sys
import numpy as np
import re

import collections
from astropy.io import ascii,fits
from astropy.table import Table, Column

"""

-------------------------------------------------------------
global variables
-------------------------------------------------------------

"""


"""

-------------------------------------------------------------
main
-------------------------------------------------------------

"""

def main(args):
    """ Main function
    """

    filter_names = args.filter_names.split(',')
    filter_depths_AB = dict(
        [(n,float(s)) for (n,s) in zip(
            filter_names, args.filter_depths_AB.split(','))])

    filter_skies_AB = dict(
        [(n,float(s)) for (n,s) in zip(
            filter_names, args.filter_skies_AB.split(','))])

    with fits.open(args.input) as tbhdu:
        data = tbhdu[1].data

    data_noise = add_noise(
        data, filter_names, filter_depths_AB, filter_skies_AB, args.r_ref)

    #if args.repeat > 1:
    #    t = Table()
    #    for c in data_noise.colnames:
    #        t[c] =  Column(np.repeat(emulated[c], args.repeat))
    #    emulated = t

    cols = []
    for d in data_noise:
        cols.append(
            fits.Column(name=d, array=data_noise[d], format='E'))

    t = fits.BinTableHDU.from_columns(
        data.columns+fits.ColDefs(cols))
    t.writeto(args.output, overwrite=True)

    return

"""

-------------------------------------------------------------
Main functions
-------------------------------------------------------------


"""




"""

-------------------------------------------------------------
Utils
-------------------------------------------------------------


"""


# function definitions
def AB_to_muJy(AB):
    return pow(10.0, -0.4*(AB - 23.9))

def muJy_to_AB(muJy):
    return -2.5*np.log10(muJy)+23.9

def compute_err(f, f_err_bgk_ref, r, r_ref, sky_AB):
    """ Compute the total flux error from the
    sky background and flux photon count.

    INPUT:
    - f_err_bgk_ref: the background flux error
    for typical Euclid galaxies
    - r_ref: is the typical size of sources
    - sky_AB is the AB magnitude of the
    background per arcsec^{-2}

    OUPUT:
    - flux error
    """

    # background flux error adjusted to
    # the object aperture
    f_err_bgk = f_err_bgk_ref * r / r_ref

    # background flux adjusted to
    # the object aperture
    f_bgk = AB_to_muJy(sky_AB)*np.pi*r**2

    # object photon error
    f_err_obj = f_err_bgk*np.sqrt(f/f_bgk)

    # total error
    return np.sqrt(f_err_bgk**2+f_err_obj**2)

def compute_mags(f, f_err, flux_limit, f_true):
    """ Compute mags and mag errors from flux
    and flux errors.

    When the flux is negative or close to
    1-sigma detection, put -1.
    """

    N = len(f)
    mag = np.zeros(N)-1.0
    magErr = np.zeros(N)-1.0

    # detected = (f > fluxLim/50.0) & (f > 0.0) & (f_true > 0.0)
    detected = (f > 0.0) & (f_true > 0.0)
    mag[detected] = -2.5*np.log10(f[detected])+23.9
    magErr[detected] = 1.0857*f_err[detected]/f_true[detected]
    # mag = muJy_to_AB(f)
    # magErr = 1.0857*f_err/f

    return mag, magErr

def add_noise(s, filters, depths_AB, sky_AB, r_ref):
    """ Add simulated fluxes to data.
    """

    result = {}
    for f in filters:

        # compute unique flux error per band from background rms
        fluxLim = pow(10.0, (-depths_AB[f]+23.9)/2.5)

        result[f+'_obs_err'] = compute_err(
            s[f], fluxLim/10.0, s['radius'],
            r_ref, sky_AB[f])
        result[f+'_obs'] = s[f] + np.random.normal(0.0, result[f+'_obs_err'])
        result[f+'_obs_mag'], result[f+'_obs_mag_err'] = compute_mags(
            result[f+'_obs'], result[f+'_obs_err'], fluxLim, s[f])

    return result

"""

-------------------------------------------------------------
Main
-------------------------------------------------------------


"""

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    # parser.add_argument('option', help="Which quantity to plot")
    parser.add_argument('-i', '--input',  default=None, help='input file')
    parser.add_argument('-o', '--output', default=None, help='output file')
    parser.add_argument('-seed', default=None, type=int, help='random seed')
    parser.add_argument(
        '-r_ref', default=None, type=float,
        help='reference source size on the sky in which depths are defined')

    parser.add_argument(
        '-filter_names', default=None, type=str, help='filter names')
    parser.add_argument(
        '-filter_skies_AB', default=None, type=str,
        help='sky brightness (AB) in filters')
    parser.add_argument(
        '-filter_depths_AB', default=None, type=str, help='filter depths (AB)')
    parser.add_argument(
        '-repeat', default=1, type=long, help='Number of repetitions')

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(seed = args.seed)




    main(args)
