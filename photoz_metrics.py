#!/usr/bin/env python

"""
Jean coupon - 2018
scripts to analyse photo-z results
"""

from __future__ import print_function, division

import os, sys
import numpy as np
import re
import collections

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import ascii,fits
from astropy.table import Table, Column

from matplotlib import rc

"""

-------------------------------------------------------------
global variables
-------------------------------------------------------------

"""

FIGX = 12.0
FIGY = 8.0
BINS_TOMO = np.array([
    0.2, 0.45, 0.55, 0.70, 0.80,
    0.9, 1.0, 1.15, 1.35, 1.65, 2.0
    ])


"""

-------------------------------------------------------------
main
-------------------------------------------------------------

"""

def main(args):
    """ Main function
    """

    rc('font',**{'size': 25, 'family':'serif', 'serif':['Palatino']})
    rc('text', usetex=True)

    # read data
    data = read_data(
        args.input, input_type=args.input_type, select=args.select)

    # loop over tasks
    print_message('Plotting...')
    if args.task == 'scatter':
        rc('figure', figsize = (FIGY, FIGY))
        plot_scatter(
            data['z_ref'], data['z'],
            args.zmin, args.zmax,
            file_output=args.output,
            title=args.title,
            stat_file_output=args.stats_output,
            density=args.density)

    elif args.task == 'PDF':
        rc('figure', figsize = (FIGX, FIGY))

        z_bins = [float(s) for s in args.z_bins.split(",")]

        plot_PDF(
            data['z_ref'], data['z'],
            data['PDF'], data['PDF_bins'],
            z_bins=z_bins,
            file_output=args.output,
            title=args.title,
            stat_file_output=args.stats_output)

    else:
        raise ValueError(
            'main: task {} is not recognised.'.format(args.task))

    print_message('Done\n')

    return

"""

-------------------------------------------------------------
Main functions
-------------------------------------------------------------


"""

def plot_PDF(
        z_ref, z, PDF, PDF_bins, z_bins=BINS_TOMO, ax=None,
        file_output=None, title = '', stat_file_output = None):
    """ Plot PDFs and compute statistics
    """

    # first make sure whether a filename
    # or an axis object is given
    if ax is None and file_output is None:
        raise ValueError(
            'plot_PDF: please provide an output filename or an axis object')

    # bin centers
    # z_bins_center = (z_bins[1:]+z_bins[:-1])/2.0
    z_bins_N = len(z_bins)-1

    # figure
    if file_output is not None:
        pp = PdfPages(file_output)

    # stats
    zmin = np.zeros(z_bins_N)
    zmax = np.zeros(z_bins_N)
    N = np.zeros(z_bins_N)
    mean = np.zeros(z_bins_N)
    mean_ref = np.zeros(z_bins_N)
    f_05 = np.zeros(z_bins_N)
    f_15 = np.zeros(z_bins_N)

    # loop over redshift bins
    for i, b in enumerate(range(z_bins_N)):

        if file_output is not None:
            fig, ax = plt.subplots(1,2)

        # select redshifts
        zmin[i] = z_bins[b]
        zmax[i] = z_bins[b+1]
        s = (z_bins[b] < z) & (z < z_bins[b+1])

        # plot n(z)
        N[i], mean[i], mean_ref[i] = plot_nz(
            PDF_bins, PDF[s], z_ref[s], info=True,
            color='b', xmax=4.0, ax = ax[0])

        # plot sum of PDFs
        f_05[i], f_15[i] = plot_PDF_dist(
            PDF_bins, PDF[s], z[s], z_ref[s],
            np.mean(z[s]), xlim=(-0.5,+0.5), ax = ax[1])

        # save figure
        if file_output is not None:
            # fig.set_tight_layout(True)
            fig.subplots_adjust(wspace=0.4)
            fig.suptitle(title+'${0:.2f}<z<{1:.2f}$'.format(zmin[i], zmax[i]))
            pp.savefig()

    if stat_file_output is not None:
        stats = collections.OrderedDict()
        stats['zmin'] = zmin
        stats['zmax'] = zmax
        stats['N'] = N
        stats['zmean'] = mean
        stats['zmean_ref'] = mean_ref
        stats['f_05'] = f_05
        stats['f_15'] = f_15

        # Table(stats).write(stat_file_output, overwrite=True, names=stats.keys())
        write_dict(stat_file_output, stats)

    if file_output is not None:
        pp.close()

    return


def plot_scatter(
        z_ref, z, zmin, zmax, ax=None,
        file_output=None, title = '',
        stat_file_output = None, density = False):
    """ Plot scatter and compute statistics
    """

    # first make sure whether a filename
    # or an axis object is given
    if ax is None and file_output is None:
        raise ValueError(
            'plot_scatter: please provide an output filename or an axis object')

    if file_output is not None:
        fig, ax = plt.subplots(1,1)

    ax.set_xlim([zmin, zmax])
    ax.set_ylim([zmin, zmax])

    # scatter plot
    if density:
        hh, locx, locy = np.histogram2d(
            z_ref, z, range=[[zmin, zmax],[zmin, zmax]], bins=[150, 150])
        hh[hh == 0.0] = np.nan
        ax.imshow(
            np.log(hh.T), origin="lower", cmap='Blues',
            extent=np.array([[zmin, zmax],[zmin, zmax]]).flatten(),
            aspect='auto', interpolation="nearest")
    else:
        ax.scatter(
            z_ref, z, s=0.05,
            c='blue', label=r'', marker = 'o')

    # stats
    stats = get_stats(z, z_ref, [zmin, zmax])

    # write stats in file
    if stat_file_output is not None:
        Table(stats).write(stat_file_output, overwrite=True, names=stats.keys())

    stats_string  = '$N_\mathrm{{gals}} = {0}$'.format(stats['N'][0])
    stats_string += '\n$\sigma = {0:5.3f}\\times(1+z)$'.format(
        stats['sigma'][0])
    stats_string += '\n$\eta = {0:5.2f}\%$'.format(stats['eta'][0])
    stats_string += '\n$\mathrm{{bias}} = {0:5.3f}\\times(1+z)$'.format(
        stats['bias'][0])

    ax.text(
        0.03*(zmax-zmin)+zmin, 0.73*(zmax-zmin)+zmin,
        stats_string, fontsize = 'x-small')

    # red lines
    x = np.arange(zmin, zmax, 0.01)
    ax.plot(x, x, 'r-', lw=1)
    ax.plot(x, x + 0.15*(1+x), 'r:', lw=1)
    ax.plot(x, x - 0.15*(1+x), 'r:', lw=1)

    # axis labels
    ax.set_xlabel(r'$z_\mathrm{ref}$')
    ax.set_ylabel(r'$z_\mathrm{phot}$')

    # title
    ax.set_title(title)

    # save figure
    if file_output is not None:
        # fig.set_tight_layout(True)
        fig.savefig(file_output)

    return

"""

-------------------------------------------------------------
Utils
-------------------------------------------------------------


"""


def plot_PDF_dist(PDFBins, PDF, zp, zs, z0, ax=None, xlim=(-2.0,2.0)):
    """ plot the summed PDF distribution
    """

    # Compute the PDF(z-z_true) sum
    dist, bins = PDF_dist(PDFBins, PDF, -4.0, +4.0, zs)

    # mode, mean and standard deviation of the distribution
    N = len(PDF)
    mode = bins[np.argmax(dist)]
    mean = trapz_boundaries(bins, bins*dist, -4.0, 4.0)/(1.0+z0)
    std_dev = np.sqrt(
        trapz_boundaries(bins, np.square(bins-mean)*dist, -4.0, 4.0))/(1.0+z0)

    # fraction of probability within 0.05(1+z) and 0.15(1+z)
    z_68 = 0.05*(1.0+z0)
    z_99 = 0.15*(1.0+z0)

    f_05 = trapz_boundaries(bins-mode, dist, -z_68, +z_68)
    f_15 = trapz_boundaries(bins-mode, dist, -z_99, +z_99)

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(bins, dist, color='b', lw=3)
    ax.fill_between(bins, 0.0, dist, color='b', alpha=0.5, label='$\Sigma$',\
                     where=(-z_68+mode < bins) & (bins < +z_68+mode))
    ax.fill_between(bins, 0.0, dist, color='b', alpha=0.5, label='$n(z_{{true}}$',\
                     where= (-z_99+mode < bins) & (bins < +z_99+mode))

    xlim = ax.set_xlim(xlim); xspan = xlim[1]-xlim[0]
    ylim = ax.set_ylim(0.0,); yspan = ylim[1]-ylim[0]

    ax.set_xlabel('$\Delta z$')
    ax.set_ylabel('PDF$(z-z_{\mathrm{true}})$')

    ax.axvline(x=mean, color='b', lw=3)
    ax.axvline(x=0.0, color='r')

    info_string = '$\mathrm{{[in\,(1+z)\,unit]}}$'
    info_string += '\n$f_{{0.05}} = {0:3.2f}\%$'.format(f_05*100)
    info_string += '\n$f_{{0.15}} = {0:3.2f}\%$'.format(f_15*100)
    info_string += '\n$\langle \Delta_z \\rangle = {0:3.3f}$'.format(mean)

    ax.text(
        0.55*xspan+xlim[0], 0.70*yspan+ylim[0],
        info_string, fontsize = 'x-small')


    return f_05, f_15


def plot_nz(
        PDF_bins, PDF, z_ref, info=True, ref_nz=False,
        color='b', xmax=4.0, ax=None, xlim=(0.0,2.0)):
    """ Plots the n(z) """

    if ax is None:
        fig, ax = plt.subplots()

    # Compute the reference redshift histogram
    if ref_nz:
        hist, hist_bins = np.histogram(z_ref, bins=50, density=True)
        ax.fill_between(0.5*(hist_bins[1:]+hist_bins[:-1]), 0.0, hist,
                        color=color, alpha=0.5, label='$n(z_\mathrm{zspec})$')
        ax.plot(0.5*(hist_bins[1:]+hist_bins[:-1]), hist, color=color, lw=2)
    # Compute the PDF(z) sum
    else:
        nz, nz_bins = PDF_dist(PDF_bins, PDF, +0.0, +4.0)
        ax.fill_between(nz_bins, 0.0, nz, color=color,
                        alpha=0.5, label='$\Sigma \mathrm{PDF}$')
        ax.plot(nz_bins, nz, color=color, lw=2)

    # plot true n(z) plus info
    if info:
        mean_ref = np.mean(z_ref)
        hist, hist_bins = np.histogram(z_ref, bins=50, density=True)
        ax.fill_between(0.5*(hist_bins[1:]+hist_bins[:-1]), 0.0, hist, color='r',
                        alpha=0.5, label='$n(z_\mathrm{zspec})$')
        ax.plot(0.5*(hist_bins[1:]+hist_bins[:-1]), hist, color='r', lw=2)
        N = len(PDF)
        mean = trapz_boundaries(nz_bins, nz_bins*nz, -4.0, 4.0)
        ylim = ax.set_ylim(0.0, ); yspan = ylim[1]-ylim[0]
        xspan = xlim[1]-xlim[0]

        info_string = '$N_\mathrm{{gals}} = {0}$'.format(N)
        info_string += '\n$\langle z_\mathrm{{PDF}}\\rangle = {0:3.3f}$'.format(mean)
        info_string += '\n$\langle z_\mathrm{{spec}}\\rangle = {0:3.3f}$'.format(np.mean(z_ref))

        ax.text(0.55*xspan+xlim[0], 0.70*yspan+ylim[0], info_string, fontsize = 'x-small')

    ax.set_xlabel('$z$')
    ax.set_ylabel('$n(z)$')

    ax.legend(loc='upper right', fontsize='x-small')
    ax.set_xlim(xlim)
    if info:
        return N, mean, mean_ref
    else:
        return N, mean



def PDF_dist(x, PDF, a, b, x0=None):
    """ return the sum of normalised PDF(x-x0)
    PDF should have linearly spaced values
    """

    step = 1.e-3

    N = len(PDF)
    int_range = np.arange(a, b, step)
    result = np.zeros(len(int_range))

    # loop over objects
    if x0 is None:
        for i in range(N):
            if np.sum(PDF[i]) > 0.0:
                result += np.interp(
                    int_range, x, PDF[i]/np.sum(PDF[i])/(x[1] - x[0]))
    else:
        for i in range(N):
            if np.sum(PDF[i]) > 0.0:
                result += np.interp(
                    int_range, x-x0[i], PDF[i]/np.sum(PDF[i])/(x[1] - x[0]))

    # return normalised sum and bins
    if N > 1:
        return result/N, int_range
    else:
        return 0.0, int_range


def get_stats(x, y, bins):
    """ Compute statistics of x vs y.
    Return (outlier free) sigma, outlier rate
    and bias
    """

    bins_N = len(bins)-1

    N = np.zeros(bins_N)
    sigma = np.zeros(bins_N)
    eta = np.zeros(bins_N)
    bias = np.zeros(bins_N)

    for i in range(0, bins_N):
        select = (bins[i] < x) & (x < bins[i+1])
        dist = (x[select]-y[select])/(1+y[select])
        n_outliers = np.sum(abs(dist) > 0.15)

        N[i] = len(dist)
        if N[i] > 0:
            sigma[i] = 1.48*np.median(np.abs(dist))
            eta[i] = n_outliers/N[i] * 100.0
            bias[i] = np.median(dist)

    result = collections.OrderedDict()
    result['zmin'] = bins[:-1]
    result['zmax'] = bins[1:]
    result['N'] = N
    result['sigma'] = sigma
    result['eta'] = eta
    result['bias'] = bias

    return result


def read_data(
    file_input, input_type='ECLD_PHZ', select=None):
    """ Read data and return data dictionary.
    Input types:
        - ECLD_PHZ: Euclid PHZ data model
        - ... ADD YOUR DATA MODEL HERE ...

    TODO:
        - implement weighting
        - filter out NaNs

    """

    print_message('Reading {}...'.format(file_input))
    if input_type == 'ECLD_PHZ':

        # column name mapping
        col_names = {
            'z' : 'TrueRedshiftPDZ_50',
            'z_ref' : 'z_true'
        }

        # read data into astropy table format
        data = Table.read(file_input, hdu=1)

        # read PDFs
        PDF = data['TrueRedshiftPDZ'].data

        # eliminate non-finite estimates
        finite = np.isfinite(data['TrueRedshiftPDZ_50'])

        data = data[finite]
        PDF = PDF[finite]

        # PDF bins
        PDF_bins = Table.read(file_input, hdu=2)['Redshift'].data

    else:
        raise ValueError(
            'read_data: input type {} is not recognised.'.format(
                input_type))
    print_message('Done\n')

    # apply selection
    if select is not None:
        cmd = ''
        for s in re.split('(\W+)', select):
            if s in data.columns:
                cmd += 'data[\'{}\']'.format(s)
            else:
                cmd += s

        select_array = eval(cmd)

        print_message('Applying selection: '+cmd+'\n')

        data = data[select_array]
        PDF = PDF[select_array]

    result = {}
    for c in col_names:
        result[c] = data[col_names[c]]

    result['PDF'] = PDF
    result['PDF_bins'] = PDF_bins

    return result

def trapz_boundaries(x, y, a, b):
    """ return the integrated value between a and b
    use the trapeze rule """

    int_range = np.linspace(a, b, num=1000)
    f = np.interp(int_range, x, y)
    return np.sum(np.diff(int_range) * (f[:-1]+f[1:])/2.0)


def print_message(text):
    """ Print a message to the standard output_PDF
    """

    # prompt messages
    verbose = True

    if verbose:
        sys.stdout.write(text)
        sys.stdout.flush()

    return


def write_dict(file_out, dict_in):
    """Wite a dictionary into csv file """

    N_rows = len(dict_in.values()[0])
    for v in dict_in.values():
        if N_rows != len(v):
            raise ValueError(
                'write_dict: columns have not the same length.')

    with open(file_out, 'w') as f:
        f.write((','.join(dict_in.keys()))+'\n')
        for row in np.array(dict_in.values()).T:
            f.write(','.join([str(r) for r in row])+'\n')

    return


"""

-------------------------------------------------------------
Main
-------------------------------------------------------------


"""

if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser()
    parser.add_argument(
        'task', help='Metric to compute and plot [scatter, PDF]')
    parser.add_argument('-i', '--input',  default=None, help='input file')
    parser.add_argument(
        '-o', '--output', default='graph.pdf', help='output file')
    parser.add_argument(
        '-input_type',  default='ECLD_PHZ',
        help='input file type (default: ECLD_PHZ)')

    parser.add_argument(
        '-zmin', default=0.0, type=float,
        help='minimum redshift for the scatter plot (default: 0.0)')

    parser.add_argument(
        '-zmax', default=6.0, type=float,
        help='maximum redshift for the scatter plot (default: 6.0)')

    parser.add_argument(
        '-title', default='',
        help='title (default: None)')

    parser.add_argument(
        '-select', default=None,
        help='selection string (default: None)')

    parser.add_argument(
        '-density', action='store_true',
        help='display scatter plot as density')

    parser.add_argument(
        '-stats_output', default=None, help='stats output file')

    z_bins_string = ",".join(map(str, BINS_TOMO))
    parser.add_argument(
        '-z_bins', default=z_bins_string,
        help='redshift bins for PDF analysis (default: {})'.format(z_bins_string))

    args = parser.parse_args()

    main(args)
