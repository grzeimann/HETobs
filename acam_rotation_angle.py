#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:42:57 2024

@author: grz85
"""

import datetime
import logging
import tarfile

import argparse as ap
import matplotlib.pyplot as plt
import numpy as np
import os.path as op

from astropy.io import fits
from astropy.stats import biweight_location, biweight_midvariance, SigmaClip
from astropy.stats import biweight_location as biweight
from astropy.visualization import AsinhStretch, PercentileInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D, BiweightLocationBackground
from photutils.detection import DAOStarFinder



def setup_logging():
    '''Set up a logger for analysis with a name ``hetobs``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    log = logging.getLogger('hetobs')
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        level = logging.INFO

        fmt = logging.Formatter(fmt)

        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)

        log = logging.getLogger('hetobs')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log

def get_lrs2_filename(tarfolder):
    '''
    

    Parameters
    ----------
    tarfolder : TYPE
        DESCRIPTION.

    Returns
    -------
    name : TYPE
        DESCRIPTION.

    '''
    T = tarfile.open(tarfolder, 'r')
    flag = True
    while flag:
        a = T.next()
        try:
            name = a.name
        except:
            flag = False
        if name[-5:] == '.fits':
            return name

def get_datetime_from_filename(filename):
    ''' Rip out the date and time from HET filenaming scheme '''
    date, time = op.basename(filename).split('_')[0].split('T')
    year, month, day = [int(x) for x in (date[:4], date[4:6], date[6:])]
    hour, minute, second = [int(x) for x in (time[:2], time[2:4], time[4:6])]
    dt = datetime.datetime(year, month, day, hour, minute, second)
    return dt


def get_imager_list_by_obs(dt, tarfolder, lag=60.):
    '''
    Get List of acm, gc1, or gc2 images for observation

    Input
    -----
    dt : datetime.datetime() object
        Using the lag in seconds, return a list of filenames within the lag
        time before the observation.  In other words, all files taken
        before the observation within "lag" seconds.

    Output
    ------
    fkeep : list
        list of filenames
    '''
    T = tarfile.open(tarfolder, 'r')
    names = T.getnames()
    tdiff = np.zeros((len(names),))
    for j, name in enumerate(names):
        adt = get_datetime_from_filename(name)
        # Positive differences mean before observation
        diff = dt - adt
        tdiff[j] = diff.days * 86400. + diff.seconds
    inds = np.where((tdiff > 0.) * (tdiff < lag))[0]
    fkeep, td = ([], [])
    for ind in inds:
        name = names[ind]
        fkeep.append(fits.open(T.extractfile(name)))
        td.append(tdiff[ind])
    return fkeep, td


def reduce_acam(image):
    ''' Really simple reduction '''
    # Get overscan but skip first column because it looks off
    overscan = biweight_location(image[:, 775:])
    # Trim all of the overscan and subtracted the average value
    image = image[:, 4:] - overscan
    columns = biweight(image, axis=0)
    # Subtract the currently global variable "acam_bias" (master bias)
    image = image - columns[np.newaxis, :]
    return image[:, :775]

def measure_image_background(image):
    ''' This does not have to be sophisticated '''
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = BiweightLocationBackground()
    bkg = Background2D(image, (25, 25), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    return image-bkg.background, bkg

def find_matches(sources, xc, yc):
    ''' Matching sources using closest neighbor and clipping '''
    dx = (sources['xcentroid'][:, np.newaxis] - xc)
    dy = (sources['ycentroid'][:, np.newaxis] - yc)
    dd = np.sqrt(dx**2 + dy**2)
    ind = np.argmin(dd, axis=1)
    dxk = dx[np.arange(dx.shape[0]), ind]
    dyk = dy[np.arange(dy.shape[0]), ind]
    bvx = biweight_midvariance(dxk)
    bvy = biweight_midvariance(dyk)
    mvx = biweight_location(dxk)
    mvy = biweight_location(dyk)
    sel = np.where((np.abs(dxk-mvx) < 3.*bvx) * (np.abs(dyk-mvy) < 3.*bvy))[0]
    return sel, ind[sel]

def detect_in_image(image_sub, bkg, thresh=3., fwhm=1.8, scale=0.2709):
    ''' Key parameters here are thresh, threshold above backgroundm and fwhm'''
    fwhm_i = fwhm / scale
    d = DAOStarFinder(fwhm=fwhm_i,
                      threshold=thresh*bkg.background_rms_median,
                      exclude_border=True)
    sources = d(image_sub)
    return sources

def make_image_subplot(fig, image, wcs, title=None,
                       vval=None, use_norm=True, cmap='Greys',
                       use_projection=True):
    ''' Plotting tool for checking astrometric solution of ACAM '''
    figkwargs = {}
    kwargs = {}
    if use_projection:
        figkwargs['projection'] = wcs
    else:
        extent = [-image.shape[1]/2.*wcs.wcs.cdelt[0]*3600.,
                  image.shape[1]/2.*wcs.wcs.cdelt[0]*3600.,
                  -image.shape[0]/2.*wcs.wcs.cdelt[1]*3600.,
                  image.shape[0]/2.*wcs.wcs.cdelt[1]*3600.]
        kwargs['extent'] = extent
    fig.add_subplot(111, **figkwargs)
    kwargs['origin'] = 'lower'
    kwargs['interpolation'] = 'none'
    kwargs['cmap'] = cmap
    if use_norm:
        kwargs['norm'] = ImageNormalize(stretch=AsinhStretch())
    if vval is not None:
        kwargs['vmin'] = vval[0]
        kwargs['vmax'] = vval[1]
    else:
        vl = np.percentile(image, 3)
        vh = np.percentile(image, 99)
        ran = vh - vl
        vl = vl - 0.2*ran
        vh = vh + 1.5*ran
        kwargs['vmin'] = vl
        kwargs['vmax'] = vh
    plt.imshow(image, **kwargs)
    if title is not None:
        plt.title(title)
    ax = plt.gca()
    plt.xlim(ax.get_xlim())
    plt.ylim(ax.get_ylim())
    plt.xlabel('RA')
    plt.ylabel('Dec')
    
# HARD CODED FOR TESTING
acam_pix_scale = 0.2709

# =============================================================================
# Argument parser so the code can be run command line
# Example: python acam_rotation_angle.py 20240925 10
# =============================================================================


parser = ap.ArgumentParser(add_help=True)

parser.add_argument("date",
                    help='''Date for reduction''',
                    type=str, default=None)

parser.add_argument("observation",
                    help='''Observation ID''',
                    type=int, default=None)

# Set argv = None, when you want to run this command line and not in debug mode
argv = ['20240928', '10']
args = parser.parse_args(args=argv)

# Set the path to correct location for the LRS2 and ACM folders
rootdir = '/work/03946/hetdex/maverick'
# Right now I expect the acm folder to be tarred up
tarfolder = op.join(rootdir, args.date, 'acm', 'acm.tar')
# I expect the lrs2 observation to be tarred as well
lrs2tarfolder = op.join(rootdir, args.date, 'lrs2', 'lrs2%07d.tar' % args.observation)


# =============================================================================
# Start the main process
# =============================================================================

log = setup_logging()
log.info('Getting observation info from the lrs2 tar file')
lrs2fn = get_lrs2_filename(lrs2tarfolder)
dt = get_datetime_from_filename(lrs2fn)
log.info('Getting acam image')
acm_fits_objs, mtdiffs = get_imager_list_by_obs(dt, tarfolder)
for acm_fits_obj, mtdiff in zip(acm_fits_objs, mtdiffs):
    log.info('ACAM is %0.2fs before LRS2 observation' %mtdiff)
    log.info('Reducing acam image')
    acm_image = acm_fits_obj[0].data
    acm_image = reduce_acam(acm_image)
    acm_sub, bkg = measure_image_background(acm_image)
    sources = detect_in_image(acm_sub, bkg, thresh=5., fwhm=1.8, scale=0.2709)
    log.info('Plotting acam image')
    plt.figure(figsize=(15, 12))
    norm = ImageNormalize(acm_sub[30:-30,:], interval=PercentileInterval(98.5),
                          stretch=AsinhStretch())
    plt.imshow(acm_sub, origin='lower', norm=norm, cmap=plt.get_cmap('coolwarm'))
    plt.colorbar()
    if len(sources):
        plt.scatter(sources['xcentroid'], sources['ycentroid'], marker='x', color='k')
