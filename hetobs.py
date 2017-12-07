# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 07:30:47 2017

@author: gregz
"""
import astrometry
import datetime
import glob
import logging

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os.path as op

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import biweight_location, biweight_midvariance
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from catalog_search import queryTESS_IC
from numpy import cos, sin, tan
from photutils import Background2D, SigmaClip, BiweightLocationBackground
from photutils import DAOStarFinder

# HARD CODED FOR TESTING
acam_x_origin = 189.1 - 4.  # Trim the overscan (first four pixels)
acam_y_origin = 386.8
acam_pix_scale = 0.2709
# This can bomb because of relative link
acam_bias = fits.open('calibration/acm/masterbias_acm.fits')
acam_bias = acam_bias[0].data

#################################
# EDIT NEXT TWO LINES FOR TESTING
rootdir = '/Users/gregz/cure/lrs2_raw/20170721'
observation = ('/Users/gregz/cure/lrs2_raw/20170721/lrs2/lrs20000004'
               '/exp01/lrs2/20170721T045123.5_066LL_sci.fits')


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


def get_datetime_from_filename(filename):
    ''' Rip out the date and time from HET filenaming scheme '''
    date, time = op.basename(filename).split('_')[0].split('T')
    year, month, day = [int(x) for x in (date[:4], date[4:6], date[6:])]
    hour, minute, second = [int(x) for x in (time[:2], time[2:4], time[4:6])]
    dt = datetime.datetime(year, month, day, hour, minute, second)
    return dt


def get_obs_info(observation):
    ''' Get the position we think we are at as well as paralactic angle'''
    F = fits.open(observation)
    ra = F[0].header['TRAJCRA'] * 15.
    dec = F[0].header['TRAJCDEC']
    pa = F[0].header['PARANGLE']
    dt = get_datetime_from_filename(observation)
    return ra, dec, pa, dt


def get_imager_list_by_obs(dt, lag=300., kind='acm'):
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
    filenames = glob.glob(op.join(rootdir, kind, '2*.fits'))
    fkeep = []
    for fn in filenames:
        adt = get_datetime_from_filename(fn)
        # Positive differences mean before observation
        diff = dt - adt
        tdiff = diff.days * 86400. + diff.seconds
        # For the time difference we want files with tdiff > 0 and <= lag
        if tdiff > 0. and tdiff <= lag:
            fkeep.append(fn)
    return fkeep


def reduce_acam(image):
    ''' Really simple reduction '''
    # Get overscan but skip first column because it looks off
    overscan = biweight_location(image[:, 1:4])
    # Trim all of the overscan and subtracted the average value
    image = image[:, 4:] - overscan
    # Subtract the currently global variable "acam_bias" (master bias)
    image = image - acam_bias
    return image


def reduce_guider(kind):
    ''' Fill this in when necessary '''
    pass


def get_stars_for_obs(ra, dec, rad=25./60.):
    ''' Query the TESS catalog, but could be more complicated if need be '''
    return queryTESS_IC(ra, dec, rad)


def get_stars_for_acam(catalog, acam_TP):
    '''
    Use coordinates in catalog to select stars within 3' of ACAM central pixel
    '''
    ra0, dec0 = acam_TP.tp.wcs_pix2world(385.5, 385.5, 1)
    c = SkyCoord(np.array([ra0])*u.degree, np.array([dec0])*u.degree,
                 frame='fk5')
    cat = SkyCoord(ra=catalog['ra']*u.degree, dec=catalog['dec']*u.degree,
                   frame='fk5')
    sep = c.separation(cat)
    loc = np.where((sep.arcmin < 3.))[0]
    return catalog[loc]


def get_star_candidates_for_probe(catalog, kind):
    ''' Fill in when necessary '''
    pass


def set_initial_acam_astrometry(ra, dec, pa):
    ''' We need ra, dec, pa and hard coded values above '''
    return astrometry.Astrometry(ra0=ra, dec0=dec, pa=pa, x0=acam_x_origin,
                                 y0=acam_y_origin, x_scale=acam_pix_scale,
                                 y_scale=acam_pix_scale, kind='acam')


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


def measure_image_background(image):
    ''' This does not have to be sophisticated '''
    sigma_clip = SigmaClip(sigma=3., iters=3)
    bkg_estimator = BiweightLocationBackground()
    bkg = Background2D(image, (100, 100), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    return image-bkg.background, bkg


def detect_in_image(image_sub, bkg, thresh=5., fwhm=1.5, scale=0.2709):
    ''' Key parameters here are thresh, threshold above backgroundm and fwhm'''
    fwhm_i = fwhm / scale
    d = DAOStarFinder(fwhm=fwhm_i,
                      threshold=thresh*bkg.background_rms_median,
                      exclude_border=True)
    sources = d(image_sub)
    return sources


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


def find_astrometric_solution(ra, dec, x, y, ra0, dec0, x0, y0):
    ''' This linear algebra problem is constrained with >= 3 stars '''
    # convert to radians for formula
    rar, decr, ra0r, dec0r = [np.deg2rad(i) for i in [ra, dec, ra0, dec0]]
    # Since x0 and y0 might not be 0,0
    xs = x - x0
    ys = y - y0
    n = len(xs)
    # Standard coordinates ksi and eta
    ksi = (cos(decr) * sin(rar-ra0r) /
           (sin(decr)*sin(dec0r) + cos(decr)*cos(dec0r)*cos(rar-ra0r)))
    eta = ((sin(decr)*cos(dec0r) - cos(decr)*sin(dec0r)*cos(rar-ra0r)) /
           (sin(decr)*sin(dec0r) + cos(decr)*cos(dec0r)*cos(rar-ra0r)))
    # 2 * len(x) by 6 matrix
    # x y 1 0 0 0 # repeat for all x and y
    # 0 0 0 x y 1 # repeat for all x and y
    A = np.zeros((n * 2, 6))
    A[:n, 0], A[:n, 1], A[:n, 2] = (xs, ys, np.ones((n,)))
    A[n:, 3], A[n:, 4], A[n:, 5] = (xs, ys, np.ones((n,)))
    b = np.hstack([ksi, eta])  # size 2 * len(x)
    sol = np.linalg.lstsq(A, b)[0]
    return sol

# Run everything outside of function so I can debug in ipython
log = setup_logging()
log.info('Getting observation info')
ra0, dec0, pa, dt = get_obs_info(observation)
log.info('Gathering stars from catalog')
catalog = get_stars_for_obs(ra0, dec0, 5./60.)
log.info('Setting up astrometry')
acm_TP = set_initial_acam_astrometry(ra0, dec0, pa)
log.info('Getting acam image list')
acm_list = get_imager_list_by_obs(dt)
log.info('Getting acam stars from catalog')
cat = get_stars_for_acam(catalog, acm_TP)
# Getting x and y position of catalog stars in ACAM
xc, yc = acm_TP.tp.wcs_world2pix(cat['ra'], cat['dec'], 1)
for acm in acm_list:
    log.info('Beginning analysis for %s' % op.basename(acm))
    F = fits.open(acm)
    image = F[0].data
    log.info('Reducing image for %s' % op.basename(acm))
    image = reduce_acam(image)
    log.info('Measuring background for %s' % op.basename(acm))
    image, bkg = measure_image_background(image)
    log.info('Detecting sources for %s' % op.basename(acm))
    sources = detect_in_image(image, bkg)
    log.info('Making plot for %s' % op.basename(acm))
    fig = plt.figure(figsize=(8, 8))
    make_image_subplot(fig, image, acm_TP.tp)
    plt.scatter(sources['xcentroid'], sources['ycentroid'],
                color='g', marker='o', s=100, facecolor='none', lw=2)
    log.info('Getting offset for %s' % op.basename(acm))
    sel_src, sel_cat = find_matches(sources, xc, yc)
    x, y, ra, dec = (sources['xcentroid'][sel_src],
                     sources['ycentroid'][sel_src], cat['ra'][sel_cat],
                     cat['dec'][sel_cat])
    sol = find_astrometric_solution(ra, dec, x, y, ra0, dec0, acam_x_origin,
                                    acam_y_origin)
    log.info('Updating projection for %s' % op.basename(acm))
    # acm_TP.update_projection()
    # xc, yc = acm_TP.tp.wcs_world2pix(cat['ra'], cat['dec'], 1)
    # plt.scatter(xc, yc, marker='x', color='r')
    plt.show()
