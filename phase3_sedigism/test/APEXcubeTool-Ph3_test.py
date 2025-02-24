#!/usr/bin/env python
"""
APEX Phase 3 Tool: Creation of FITS spectral cubes compliant with ESO Phase 3 requirements.
Contact: pvenegas@eso.org, tstanke@mpe.mpg.de
Creation: 2020-07-13; Updated: 2021-08-13
"""

import os
import sys
import math
import copy
import hashlib
import datetime
import pathlib
import re

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.wcs import WCS

# Colors for printing messages
RED = "\033[31m"       # For errors
GREEN = "\033[32;3m"   # For outputs
BLUE = "\033[34;1m"    # Request inputs
MAGENTA = "\033[35m"   # Header
CYAN = "\033[36m"      # Main words
ENDCOLOR = "\033[0m"

##########################################
# CONFIGURATION: Set input parameters here.
##########################################
# Set the names of the input FITS files (must be present in the current directory)
SPEC_IN = "test_reprojected.fits"         # Science cube file
RMS_IN = "test_reprojected.fits"            # Error (rms) cube file

# LMV file association (if applicable)
ASSOCIATE_LMV = False                       # Set True if you want to associate a .lmv file
LMV_FILE = "your_cube.lmv"                  # Specify the .lmv file (if ASSOCIATE_LMV is True)

# TAP file pause: set to True if you wish to pause for manual TAP file editing
TAP_PAUSE = False

# Jy/K conversion factor handling:
# If DIFFERENT_JY_FACTOR is True, the default factor from the dictionary will be overridden.
DIFFERENT_JY_FACTOR = False
NEW_JY_FACTOR = 53.0                        # New Jy/K factor (if overriding)
NEW_JY_UNCERT = 8.0                         # New Jy/K uncertainty

# Beam efficiency correction:
# If APPLY_BEAM_EFFICIENCY is True, the beam efficiency correction will be applied.
APPLY_BEAM_EFFICIENCY = False
BEAM_EFFICIENCY = 1.0                        # Value between 0 and 1

# Flux conversion:
# Set CONVERT_TO_JY to True if you want to convert from K to Jy.
CONVERT_TO_JY = False

# Programme Identification:
# If CONFIRM_PROG_CODE is False, the new programme code NEW_PROG_CODE will be used.
CONFIRM_PROG_CODE = True
NEW_PROG_CODE = ""                          # Only used if CONFIRM_PROG_CODE is False

# Bibliographic reference (Bibcode); leave empty if none.
BIBLIO_REF = "2021MNRAS.500.3064S"

# Beam size override:
# Set NEW_BEAM_SIZE to a number (in arcsec) to override the effective beam size;
# otherwise leave as None to use the value computed from the header.
NEW_BEAM_SIZE = None

##########################################
# Dictionaries for Jy/K conversion factors
##########################################
dicc_Gen = {
    "HET230": [2016, 39, 6],
    "PI230": [2016, 44, 6],
    "HET345": [2016, 53, 8],
    "NOVA660": [2016, 137, 22],
    "GARD180": [2016, 39, 6],
    "PI230": [2017, 45, 7],
    "HET230": [2017, 40, 6],
    "HET345": [2017, 51, 8],
    "NOVA660": [2017, 110, 18],
    "GARD180": [2017, 40, 6],
    "HET460": [2017, 72, 11],
    "PI230": [2018, 46, 4],
    "SEPIA660": [2018, 63, 5],
    "GARD180": [2018, 40, 6],
    "SEPIA180": [2019, 35, 3],
    "SEPIA660": [2019, 68, 6],
}

dicc = {
    "HET230": [2016, 40, 6],
    "PI230": [2017, 45, 7],
    "HET345": [2016, 53, 8],
    "SUPERCAM": [2014, 53, 8],  # Use same value as for FLASH345 and HET345 also for SUPERCAM
    "FLASH345": [2018, 53, 4],
    "HET460": [2017, 72, 11],
    "GARD180": [2017, 40, 6],
    "NOVA660": [2017, 110, 18],
    "SEPIA180": [2019, 35, 3],
    "SEPIA660": [2019, 68, 6],
}

##########################################
# Function definitions
##########################################
def tapQuery_Het(source):
    """
    Query the ESO TAP service.
    """
    from pyvo.dal import tap
    from astropy.coordinates import SkyCoord
    from astropy.units import Quantity
    from tabulate import tabulate

    ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
    tapobs = tap.TAPService(ESO_TAP_OBS)

    long = np.int32(source.split('G')[1])
    source_left = f'G{long-1}'
    source_right = f'G{long+1}'

    print("\nQuerying the ESO TAP service at %s" % ESO_TAP_OBS)
    query = (
        "SELECT dp_id, exposure, prog_id, object, dp_tech, instrument, ra, dec \n"
        "FROM dbo.raw \n"
        "WHERE dp_id LIKE 'APEXHET.%' \n"
        "AND (prog_id LIKE '092.F-9315%' OR prog_id LIKE '193.C-0584%') \n"
        f"AND (object LIKE '{source}%' OR object LIKE '{source_left}%' OR object LIKE '{source_right}%') \n"
        "AND dp_cat = 'SCIENCE'"
    )
    print("\nQuery:\n" + query)
    res = tapobs.search(query=query, maxrec=1000)
    print("\n", res.to_table(), "\n")
    print("\nA total of " + str(len(res.to_table())) + " records were found matching the provided criteria.")
    table = res.to_table()
    filename = f"{source}.tap"
    print("\nResults has been written to: " + filename)
    with open(filename, "w") as f:
        print(tabulate(table), file=f)
    return table


hdu = fits.open('test_reprojected.fits')[0]
hdr = hdu.header
source = hdr["OBJECT"]

tapQuery_Het(source)