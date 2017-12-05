# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 07:30:47 2017

@author: gregz
"""
import astrometry
import glob
import datetime

import os.path as op
import astropy.units as u
import numpy as np

from catalog_search import queryTESS_IC
from astropy.io import fits
from astropy.stats import biweight_location
from astropy.coordinates import SkyCoord


acam_x_origin = 189.1 - 4.  # Trim the overscan (first four pixels)
acam_y_origin = 386.8
acam_pix_scale = 0.2709

acam_bias = fits.open('/Users/gregz/cure/calibration/acm/masterbias_acm.fits')
acam_bias = acam_bias[0].data

rootdir = '/Users/gregz/cure/lrs2_raw'
observation = ''


def get_datetime_from_filename(filename):
    date, time = op.basename(filename).split('_')[0].split('T')
    year, month, day = (date[:4], date[4:6], date[6:])
    hour, minute, second = (time[:2], time[2:4], time[4:6])
    dt = datetime.datetime(year, month, day, hour, minute, second)
    return dt


def get_obs_info(observation):
    ''' Get Observation Info for a filename '''
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

    Returns:
    fkeep : list
        list of filenames
    '''
    filenames = glob.glob(op.join(rootdir, kind, '2*.fits'))
    fkeep = []
    for fn in filenames:
        adt = get_datetime_from_filename(fn)
        # Positive differences mean before observation
        diff = dt - adt
        tdiff = diff.days*86400. + diff.seconds
        # For the time difference we want files with tdiff > 0 and <= lag
        if tdiff > 0. and tdiff <= lag:
            fkeep.append(fn)
    return fkeep


def reduce_acam(image):
    overscan = biweight_location(image[:, 1:4])
    image = image[:, 4:] - overscan
    image = image - acam_bias
    return image


def reduce_guider(kind):
    pass


def get_stars_for_obs(ra, dec, rad=25./60.):
    return queryTESS_IC(ra, dec, rad)


def get_stars_for_acam(catalog, acam_TP, buffer):
    ra0, dec0 = acam_TP.wcs_pix2world(385.5, 385.5, 1)
    c = SkyCoord(np.array([ra0])*u.degree, np.array([dec0])*u.degree,
                 frame='fk5')
    cat = SkyCoord(ra=catalog['ra']*u.degree, dec=catalog['dec']*u.degree,
                   frame='fk5')
    sep = c.separation(cat)
    loc = np.where((sep.arcmin < 3.))[0]
    return catalog[loc]


def get_stars_for_probe(catalog, kind):
    pass


def set_initial_acam_astrometry(ra, dec, pa):
    return astrometry.Astrometry(ra0=ra, dec0=dec, pa=pa, x0=acam_x_origin,
                                 y0=acam_y_origin, x_scale=-acam_pix_scale,
                                 y_scale=acam_pix_scale, kind='acam')


def get_astrometry_acam(init_astrom, stars):
    pass


