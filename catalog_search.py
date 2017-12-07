# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 09:28:59 2017

@author: gregz
"""

import sys
import json
from urllib import pathname2url as urlencode
import httplib
from astropy.table import Table
import numpy as np


class MakeRegionFile(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []

        s.append('# Region file format: DS9 version 4.1')
        s.append('global color=green dashlist=8 3 width=1 '
                 'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '
                 'fixed=0 edit=1 move=1 delete=1 include=1 source=1')
        s.append('fk5')
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeSource(cls, f, ra, dec, rad=2):
        s = ('circle(%0.6f, %0.6f, %0.2f")' % (ra, dec, rad))
        f.write(s + "\n")
        f.flush()


def mastQuery(request):
    ''' Mast seems to respond quickly which is great '''
    server = 'mast.stsci.edu'

    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent": "python-requests/" + version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)

    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')
    # Close the https connection
    conn.close()

    return head, content


def mastJson2Table(jsonObj):
    ''' Table() allows for numpy array's and easy calls '''
    dataTable = Table()
    for col, atype in [(x['name'], x['type']) for x in jsonObj['fields']]:
        if atype == 'string':
            atype = 'str'
        if atype == 'boolean':
            atype = 'bool'
        tcol = []
        for x in jsonObj['data']:
            col_val = x.get(col, None)
            if col_val is None:
                tcol.append(-999.)
            else:
                tcol.append(col_val)
        dataTable[col] = np.array(tcol, dtype=atype)
    return dataTable


def queryTESS_IC(ra, dec, radius):
    ''' TESS input catalog is a mash-up of 2mass, SDSS, GAIA, and more '''
    mashupRequest = {'service': 'Mast.Catalogs.Tic.Cone',
                     'params': {'ra': ra,
                                'dec': dec,
                                'radius': radius},
                     'format': 'json',
                     'pagesize': 10000,
                     'page': 1}

    headers, outString = mastQuery(mashupRequest)

    outData = json.loads(outString)

    table = mastJson2Table(outData)
    sel = np.where(table['objType'] == 'STAR')[0]
    return table[sel]
