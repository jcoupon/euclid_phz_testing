#!/usr/bin/env python

""" Jean coupon, 2016-2018

script to simulate fluxes and
flux errors

See:
http://www.cfht.hawaii.edu/Instruments/Imaging/Megacam/dietmegacam.html
http://123.physics.ucdavis.edu/week_7_files/CCD_imaging_lecture.pdf


"""

from __future__ import print_function, division

import os
import sys
import re
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column


"""

-----------------------------------------------------
Options
-----------------------------------------------------

"""

"""

-----------------------------------------------------
Data structures
-----------------------------------------------------

"""

class filterClass(object):
    """This class is a container for filter properties
    to compute flux errors.

    Pixels counts are also called
    "Data Numbers" (DN) or "Analog-to-Digital Units" (ADU)
    """

    def __init__(
            self, name='', gain=1.62, readout_noise=5.0,
            pixel_scale=0.187, k=0.150, zero_point=27.00, sky=22.0):

        """ Initialisation

        name: name of the filter
        gain: gain of the camera [e-/ADU]
        readout_noise: noise in electron produced by the amplifier [e-/pixel]
        pixel_scale: size of one pixel in arsec/pix
        k: airmass term
        zero_point = sensitivity of the camera for a given filter [e-/sec],
        expressed in mag (count per second at 1 airmass)
        """

        self.name = name
        self.gain = gain
        self.readout_noise = readout_noise
        self.pixel_scale = pixel_scale
        self.k = k
        self.zero_point = zero_point
        self.sky = sky

    def get_flux_muJy(self, fs, t_exp, N_exp=1, airmass=1.2):
        """ convert flux in electron per second to muJy

        fs: flux in electron
        t_exp: exposure time in seconds
        N_exp: number of exposures (default: 1)
        airmass: air mass between source and telescope (default: 1.2)
        """

        # flux of the source in e-/s within the aperture
        f = fs/(t_exp*N_exp*pow(10.0,
            -0.4*(23.9 - self.zero_point+self.k*(airmass-1.0))))

        return f


    def get_flux_electron(self, flux, t_exp, N_exp=1, airmass=1.2, unit='muJy'):
        """ Flux in electron per second

        flux: flux
        t_exp: exposure time in seconds
        N_exp: number of exposures (default: 1)
        airmass: air mass between source and telescope (default: 1.2)
        unit: unit of the flux (default: 'muJy')
        """

        f = to_muJy(flux, unit)

        # total flux of the source in e-/s within the aperture
        fs = t_exp*N_exp*f*pow(
            10.0, -0.4*(23.9 - self.zero_point+self.k*(airmass-1.0)))

        return fs

    def get_flux_electron_error(
            self, flux, t_exp, N_exp=1, airmass=1.2, diameter=2.0, unit='muJy'):
        """ Flux error in electron per second

        flux: flux
        t_exp: exposure time in seconds
        N_exp: number of exposures (default: 1)
        airmass: air mass between source and telescope (default: 1.2)
        diameter: circular (default: 2.0")
        unit: unit of the flux (default: 'muJy')
        """

        f = to_muJy(flux, unit)

        # total flux of the source in e-/s within the aperture
        fs = self.get_flux_electron(f, t_exp, N_exp=N_exp, airmass=airmass)

        # area in arcsec^2
        area = pow(diameter/2.0, 2.0)*np.pi

        # total flux background flux in e-/s within the aperture
        fb = t_exp*N_exp*pow(10.0, -0.4*(
            self.sky - self.zero_point + self.k*(airmass-1.0)))*area

        # read out variance within the aperture
        var_readout = N_exp*self.readout_noise**2*area/self.pixel_scale**2.0

        return np.sqrt(fs+fb+var_readout)


    def get_flux_error(
            self, flux, t_exp, N_exp=1, airmass=1.2, diameter=2.0, unit='muJy'):
        """ Return the flux error
        """

        # get the flux error in electrons
        fe_err = self.get_flux_electron_error(
                flux, t_exp, N_exp=N_exp, airmass=airmass,
                diameter=diameter, unit=unit)

        # convert back to muJy flux
        f = self.get_flux_muJy(fe_err, t_exp, N_exp=N_exp, airmass=airmass)

        return from_muJy(f, unit=unit)

    def get_SNR(
        self, flux, t_exp, N_exp=1, airmass=1.2, diameter=2, unit='muJy'):
        """ Signal-to-noise ratio

        AB: AB magnitude
        t_exp: exposure time in seconds
        N_exp: number of exposures (default: 1)
        airmass: air mass between source and telescope (default: 1.2)
        diameter: circular (default: 2.0")
        """

        flux_electron = self.get_flux_electron(
            flux, t_exp, N_exp=N_exp, airmass=airmass, unit=unit)
        flux_electron_err = self.get_flux_electron_error(
            flux, t_exp, N_exp=N_exp, airmass=airmass,
            diameter=diameter, unit=unit)

        return flux_electron/flux_electron_err


def to_muJy(flux, unit):
    """ Convert flux into micro jansky.
    """

    if unit == 'muJy':
        return flux
    else:
        raise ValueError(
            'to_muJy: unit {} for input flux is not recognised.'.format(unit))

def from_muJy(flux, unit):
    """ Convert micro Jansky flux into other kind.
    """

    if unit == 'muJy':
        return flux
    else:
        raise ValueError(
            'to_muJy: unit {} for input flux is not recognised.'.format(unit))


def AB_to_muJy(AB):
    """ Convert AB magnitudes
    into fluxes.

    args:
        AB magnitudes

    returns:
        muJy: flux in micro jansky
    """

    finite = np.isfinite(AB)

    if isinstance(AB, (list, tuple, np.ndarray)):
        # only convert finite values
        result = np.full(np.shape(AB), np.nan)
        result[finite] = pow(10.0, -0.4*(AB[finite] - 23.9))
    else:
        result = np.nan
        if finite:
            result = pow(10.0, -0.4*(AB - 23.9))

    return result

def muJy_to_AB(muJy):
    """ Convert fluxes into AB
    magnitudes

    args:
        muJy: flux in micro jansky

    returns:
        AB magnitudes
    """

    # only convert finite values
    finite = np.isfinite(muJy)

    # only convert non zero fluxes
    finite[finite] &= muJy[finite]>1.e-10

    if isinstance(muJy, (list, tuple, np.ndarray)):

        result = np.full(np.shape(muJy), np.nan)

        # AB mags
        result[finite] = -2.5*np.log10(muJy[finite])+23.9
    else:
        result = np.nan
        if finite:
            result = -2.5*np.log10(muJy)+23.9

    return result

"""

-----------------------------------------------------
Filter properties
-----------------------------------------------------

"""


"""
CFHT MegaCam http://cfht.hawaii.edu/Instruments/Imaging/MegaCam/generalinformation.html
TODO


CFHT WirCam http://www.cfht.hawaii.edu/Instruments/Imaging/WIRCam/dietWIRCam.html
TODO

Subaru HSC https://www.subarutelescope.org/Observing/Instruments/HSC/index.html
sky background computed from ETC https://hscq.naoj.hawaii.edu/cgi-bin/HSC_ETC/hsc_etc.cgi

g: dark: 22.05, grey:21.29 bright:20.13
r2: dark: 21.21, grey:20.99 bright:20.25
i2: dark: 20.20, grey:20.09 bright:19.84
z: dark: 19.88, grey:19.70 bright:19.31
Y: dark: 18.72, grey:18.55 bright:18.19

"""

FILTERS={

    # Euclid - From Jean-Charles and Jerome Amiaux see https://mail.google.com/mail/u/0/?ui=2&shva=1#inbox/1603671070901d92
    # TODO: add telescope background for Euclid
    'Euclid_VIS': filterClass(
        name='Euclid_VIS', gain=3.1, readout_noise=4.5, pixel_scale=0.1,
        k=0.0, zero_point=25.51, sky=22.33),
    'Euclid_Y': filterClass(
        name='Euclid_Y', gain=1.3, readout_noise=7.7, pixel_scale=0.3,
        k=0.0, zero_point=24.65, sky=22.10),
    'Euclid_J': filterClass(
        name='Euclid_J', gain=1.3, readout_noise=7.7, pixel_scale=0.3,
        k=0.0, zero_point=24.73, sky=22.11),
    'Euclid_H': filterClass(name='Euclid_H', gain=1.3, readout_noise=7.7,
        pixel_scale=0.3, k=0.0, zero_point=24.73, sky=22.28),

    # CFHT MegaCam - old filters ("S" = short)
    'MegaCam_uS': filterClass(
        name='MegaCam_uS', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.350, zero_point=25.74, sky=22.70), # dark
    'MegaCam_gS': filterClass(
        name='MegaCam_gS', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.150, zero_point=27.00, sky=22.0),
    'MegaCam_rS': filterClass(
        name='MegaCam_rS', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.100, zero_point=26.50, sky=21.30), # dark
    'MegaCam_iS': filterClass(
        name='MegaCam_iS', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.040, zero_point=26.38, sky=20.30), # grey
    'MegaCam_zS': filterClass(
        name='MegaCam_zS', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.030, zero_point=25.34, sky=19.40), # grey

    # CFHT MegaCam - new filters
    'MegaCam_u': filterClass(
        name='MegaCam_u', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.350, zero_point=25.78, sky=22.70), # dark
    'MegaCam_g': filterClass(
        name='MegaCam_g', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.150, zero_point=27.11, sky=22.00), # dark
    'MegaCam_r': filterClass(
        name='MegaCam_r', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.100, zero_point=26.74, sky=21.30), # dark
    'MegaCam_i': filterClass(
        name='MegaCam_i', gain=1.62, readout_noise=5.0, pixel_scale=0.187,
        k=0.040, zero_point=26.22, sky=20.30), # grey
    'MegaCam_z': filterClass(
        name='MegaCam_z', gain=1.62, readout_noise=5.0,
        pixel_scale=0.187, k=0.030, zero_point=25.02, sky=19.40), # grey

    # CFHT Wircam
    'WirCam_Ks': filterClass(
        name='WirCam_Ks', gain=3.7, readout_noise=30.0, pixel_scale=0.307,
        k=0.05, zero_point=26.371, sky=15.731),

    # Subaru HSC
    'HSC_g': filterClass(
        name='HSC_g', gain=3.0, readout_noise=4.5, pixel_scale=0.17,
        k=0.0, zero_point=29.00, sky=22.05), # dark
    'HSC_r2': filterClass(
        name='HSC_r2', gain=3.0, readout_noise=4.5, pixel_scale=0.17,
        k=0.0, zero_point=29.10, sky=21.21), # dark
    'HSC_i2': filterClass(
        name='HSC_i2', gain=3.0, readout_noise=4.5, pixel_scale=0.17,
        k=0.0, zero_point=28.70, sky=20.20), # dark
    'HSC_z': filterClass(
        name='HSC_z', gain=3.0, readout_noise=4.5, pixel_scale=0.17,
        k=0.0, zero_point=27.70, sky=19.70), # grey
    'HSC_Y': filterClass(
        name='HSC_Y', gain=3.0, readout_noise=4.5, pixel_scale=0.17,
        k=0.0, zero_point=27.40, sky=18.55) # grey
}


"""

-----------------------------------------------------
main
-----------------------------------------------------

"""


def main(args):
    function = getattr(sys.modules[__name__], args.option)(args)
    return

"""

-----------------------------------------------------
Main functions
-----------------------------------------------------

"""


def test(args):
    """
    This routines tests the flux error
    value compared with the online ETCs

    ETCs:
    CFHT-MegaCam: http://etc.cfht.hawaii.edu/mp/
    CFHT-Wircam: http://etc.cfht.hawaii.edu/wc/
    Subaru-HSC: https://hscq.naoj.hawaii.edu/cgi-bin/HSC_ETC/hsc_etc.cgi

    for CFHT, select "galaxy" with fixed aperture
    with 1" radius (=2" diameter) and select seeing = 0.01
    to make sure 100% of the flux is within the aperture

    """

    AB = 20.0
    t_exp = 100
    N_exp = 20

    print('MegaCam_gS, SNR=', FILTERS['MegaCam_gS'].get_SNR(
        AB_to_muJy(AB), t_exp, N_exp=N_exp))

    AB = 20.0
    t_exp = 2000
    N_exp = 1

    print('WirCam_Ks, SNR=', FILTERS['WirCam_Ks'].get_SNR(
        AB_to_muJy(AB), t_exp, N_exp=N_exp))

    AB = 20.0
    t_exp = 300
    N_exp = 1

    print('HSC_g, SNR=', FILTERS['HSC_g'].get_SNR(
        AB_to_muJy(AB), t_exp, N_exp=N_exp))

    AB = 24.5
    t_exp = 565
    N_exp = 3

    print('Euclid_VIS, SNR=', FILTERS['Euclid_VIS'].get_SNR(
        AB_to_muJy(AB), t_exp, N_exp=N_exp, diameter=2.0)) #1.3)

    # AB = 24.00
    # diameter = 1.016

    AB = 23.00
    diameter = 2.0

    t_exp = 95 # J: 82, H: 105
    N_exp = 3

    print('Euclid_Y, SNR=', FILTERS['Euclid_Y'].get_SNR(
        AB_to_muJy(AB), t_exp, N_exp=N_exp, diameter=diameter))

    t_exp = 91 # J: 82, H: 105
    N_exp = 3

    print('Euclid_J, SNR=', FILTERS['Euclid_J'].get_SNR(
        AB_to_muJy(AB), t_exp, N_exp=N_exp, diameter=diameter))

    t_exp = 58 # J: 82, H: 105
    N_exp = 3

    print('Euclid_H, SNR=', FILTERS['Euclid_H'].get_SNR(
        AB_to_muJy(AB), t_exp, N_exp=N_exp, diameter=diameter))

    return

def sky_brightness_HSC(args):
    """ computes the sky brightness in
    mag/arcsec^2 from Subaru ETC https://hscq.naoj.hawaii.edu/cgi-bin/HSC_ETC/hsc_etc.cgi
    """

    fs = np.array([17.41, 35.08, 101.64])/FILTERS['HSC_g'].pixel_scale**2.0
    print('g: dark: {0:.2f}, grey:{1:.2f} bright:{2:.2f}'.format(
        *FILTERS['HSC_g'].get_AB_from_flux_electron(fs, 1.0)))

    fs = np.array([41.33, 50.48,100.14])/FILTERS['HSC_r2'].pixel_scale**2.0
    print('r2: dark: {0:.2f}, grey:{1:.2f} bright:{2:.2f}'.format(
        *FILTERS['HSC_r2'].get_AB_from_flux_electron(fs, 1.0)))

    fs = np.array([72.92, 80.41, 100.98])/FILTERS['HSC_i2'].pixel_scale**2.0
    print('i2: dark: {0:.2f}, grey:{1:.2f} bright:{2:.2f}'.format(
        *FILTERS['HSC_i2'].get_AB_from_flux_electron(fs, 1.0)))

    fs = np.array([38.87, 45.60, 65.33])/FILTERS['HSC_z'].pixel_scale**2.0
    print('z: dark: {0:.2f}, grey:{1:.2f} bright:{2:.2f}'.format(
        *FILTERS['HSC_z'].get_AB_from_flux_electron(fs, 1.0)))

    fs = np.array([85.37, 100.15, 140.06])/FILTERS['HSC_Y'].pixel_scale**2.0
    print('Y: dark: {0:.2f}, grey:{1:.2f} bright:{2:.2f}'.format(
        *FILTERS['HSC_Y'].get_AB_from_flux_electron(fs, 1.0)))

    return

def sim_HSC_Wide(args):
    """ Simulates HSC Wide photometry
    """

    verbose = True

    # options

    # filter and column names
    col_names = ['hsc_g', 'hsc_r', 'hsc_i', 'hsc_z', 'hsc_y']
    filter_names = ['HSC_g', 'HSC_r2', 'HSC_i2', 'HSC_z', 'HSC_Y']

    # exposure times
    t_exp = [10.0*60.0, 10.0*60.0, 20.0*60.0, 20.0*60.0, 20.0*60.0]

    # read data as a table
    table = Table.read( args.input, hdu=1)

    # loop over filters
    for c,f,t in zip(col_names, filter_names, t_exp):

        # flux in muJy
        flux = table[c]
        flux_error = FILTERS[f].get_flux_error(flux, t)
        SNR = flux/flux_error

        # magnitude errors
        mag_error = np.full(np.shape(flux), np.nan)
        pos = flux > 0.0
        mag_error[pos] = 1.0857/SNR[pos]

        table[f+'_obs'] = np.random.normal(flux, flux_error)
        table[f+'_obs_err'] = flux_error
        table[f+'_mag_obs'] = muJy_to_AB(table[f+'_obs'])
        table[f+'_mag_obs_err'] = mag_error

    table.write(args.output, overwrite=True)

    #fileIn = fits.open(fileInName)
    #tbhdu = fits.HDUList([fileIn[0], fits.BinTableHDU.from_columns(
    #    fileIn[1].columns + fits.ColDefs(cols))])

    #tbhdu.writeto(args.output, clobber=True)

    return

def sim_HSC_training(args):
    """ Reference sample. HSC fluxes and magnitudes
    with no errors.
    """

    verbose = True

    # options

    # filter and column names
    col_names = ['hsc_g', 'hsc_r', 'hsc_i', 'hsc_z', 'hsc_y']
    filter_names = ['HSC_g', 'HSC_r2', 'HSC_i2', 'HSC_z', 'HSC_Y']

    # exposure times
    t_exp = [10.0*60.0, 10.0*60.0, 20.0*60.0, 20.0*60.0, 20.0*60.0]

    # read data as a table
    table = Table.read( args.input, hdu=1)

    # loop over filters
    for c,f,t in zip(col_names, filter_names, t_exp):

        # flux in muJy
        flux = table[c]
        AB = muJy_to_AB(flux)

        table[f+'_obs'] = flux
        table[f+'_obs_err'] = np.zeros(len(flux))
        table[f+'_mag_obs'] = AB
        table[f+'_mag_obs_err'] = np.zeros(len(flux))

    table.write(args.output, overwrite=True)

    #fileIn = fits.open(fileInName)
    #tbhdu = fits.HDUList([fileIn[0], fits.BinTableHDU.from_columns(
    #    fileIn[1].columns + fits.ColDefs(cols))])

    #tbhdu.writeto(args.output, clobber=True)

    return

"""

-----------------------------------------------------
Main
-----------------------------------------------------

"""

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('option', help="Which action")
    parser.add_argument('-i', '--input', default=None, help='input file(s)')
    parser.add_argument('-o', '--output', default=None, help='output file(s)')
    parser.add_argument(
        '-seed', default=20091982, type=long,
        help='Random seed. Default: 20091982')

    args = parser.parse_args()

    np.random.seed(seed=args.seed)

    main(args)
