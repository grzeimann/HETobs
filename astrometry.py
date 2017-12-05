# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 09:19:46 2017

@author: gregz
"""

import logging
import sys

from pyhetdex.het.fplane import FPlane
from numpy import cos, sin, deg2rad
from astropy import wcs


class Astrometry:
    def __init__(self, ra0, dec0, pa, x0, y0, x_scale=-1., y_scale=1.,
                 sys_rot=1.3, fplane_file=None, kind='fplane'):
        self.setup_logging()
        self.ra = ra0
        self.dec = dec0
        self.dra = 0.
        self.ddec = 0.
        self.x0 = x0
        self.y0 = y0
        self.pa = pa
        self.sys_rot = sys_rot
        self.fplane_file = fplane_file
        if self.fplane_file is None:
            self.log.info('No fplane file given.')
            self.log.info('Some functions will be unavailable until'
                          'an fplane is given.')
        else:
            self.fplane = FPlane(self.fplane_file)
        self.kind = kind
        self.set_effective_rotation()

        # Building tangent plane projection with scale 1"
        self.tp = self.setup_TP(self.ra0, self.dec0, self.rot, self.x0,
                                self.y0, x_scale=x_scale, y_scale=y_scale)
        self.tp_ifuslot = None

    def setup_TP(self, ra0, dec0, rot, x0=0.0, y0=0.0, x_scale=-1.,
                 y_scale=1.):
        ARCSECPERDEG = 1.0/3600.0

        # make a WCS object with appropriate FITS headers
        tp = wcs.WCS(naxis=2)
        tp.wcs.crpix = [x0, y0]  # central "pixel"
        tp.wcs.crval = [ra0, dec0]  # tangent point
        tp.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        # pixel scale in deg.
        tp.wcs.cdelt = [ARCSECPERDEG * x_scale, ARCSECPERDEG * y_scale]

        # Deal with PA rotation by adding rotation matrix to header
        rrot = deg2rad(rot)
        # clockwise rotation matrix
        tp.wcs.pc = [[cos(rrot), sin(rrot)], [-1.0*sin(rrot), cos(rrot)]]
        return tp

    def set_polynomial_platescale(self):
        self.tp.wcs.a_0_0 = 0.177311
        self.tp.wcs.a_1_0 = -8.29099e-06
        self.tp.wcs.a_2_0 = -2.37318e-05
        self.tp.wcs.b_0_0 = 0.177311
        self.tp.wcs.b_0_1 = -8.29099e-06
        self.tp.wcs.b_0_2 = -2.37318e-05
        self.tp.wcs.a_order = 2
        self.tp.wcs.b_order = 2

    def setup_logging(self):
        '''Set up a logger for analysis with a name ``shot``.

        Use a StreamHandler to write to stdout and set the level to DEBUG if
        verbose is set from the command line
        '''
        self.log = logging.getLogger('astrometry')
        if not len(self.log.handlers):
            fmt = '[%(levelname)s - %(asctime)s] %(message)s'
            level = logging.INFO

            fmt = logging.Formatter(fmt)

            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            handler.setLevel(level)

            self.log = logging.getLogger('astrometry')
            self.log.setLevel(logging.DEBUG)
            self.log.addHandler(handler)

    def set_effective_rotation(self):
        # Making rotation from the PA
        # TODO: re-derive this formula 
        if self.kind == 'fplane':
            self.rot = 360. - (90. + self.pa + self.sys_rot)
        elif self.kind == 'acam':
            self.rot = -self.pa + self.sys_rot
        else:
            self.log.error('"kind" was not set to available options.')
            self.log.info('Available options are: %s and %s' % ('fplane',
                                                                'acam'))
            self.log.info('Next time please choose one of the options above.')
            self.log.info('Exiting due to error.')
            sys.exit(1)

    def update_projection(self):
        self.set_effective_rotation()
        # Building tangent plane projection with scale 1"
        # dra = self.dra / 3600. / np.cos(np.deg2rad(self.dec))
        # ddec = self.ddec / 3600.
        self.tp = self.setup_TP(self.ra0, self.dec0, self.rot, self.x0,
                                self.y0)

    def get_ifuslot_ra_dec(self, ifuslot):
        if self.fplane is None:
            return None
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y
        return self.tp.xy2raDec(ifu.y, ifu.x)

    def get_ifuspos_ra_dec(self, ifuslot, x, y):
        if self.fplane is None:
            return None
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y
        return self.tp.wcs_pix2world(ifu.y + x, ifu.x + y, 1)

    def get_ifuslot_projection(self, ifuslot, imscale, crx, cry):
        if self.fplane is None:
            return None
        ra, dec = self.get_ifuslot_ra_dec(ifuslot)
        self.tp_ifuslot = self.setup_TP(ra, dec, self.rot, crx, cry,
                                        x_scale=-imscale, y_scale=imscale)

    def convert_ifuslot_xy_to_new_xy(self, x, y, wcs):
        if self.tp_ifuslot is None:
            self.log.error('You have not setup the ifuslot projection yet.')
            self.log.error('To do so, call '
                           '"get_ifuslot_projection(ifuslot, imscale')
            return None
        ra, dec = self.tp_ifuslot.wcs.wcs_pix2world(x, y, 1)
        return wcs.wcs_world2pix(ra, dec, 1)
