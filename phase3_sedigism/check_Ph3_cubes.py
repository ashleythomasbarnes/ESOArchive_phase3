#!/usr/bin/env python
"""
APEX Phase 3 Tool: Creation of FITS spectral cubes compliant with ESO Phase 3 requirements.
Contact: pvenegas@eso.org, tstanke@mpe.mpg.de
Creation: 2020-07-13; Updated: 2021-08-13
"""

##########################################
# TO DO 
# - Pull this from TAP but need ISTs for this... (currently hardcoded)
##########################################

# Imports
import os # For file handling 
import sys # For system handling
import math # For mathematical operations
import copy # For copying objects
import hashlib # For hashing

import numpy as np # For numerical operations
from astropy.io import fits # For FITS file handling
from astropy.time import Time # For time handling
from astropy import units as u # For units
from astropy.table import Table # For table handling
# from astropy.wcs import WCS # For World Coordinate System
from astropy.coordinates import SkyCoord # For sky coordinates
from astropy.io.votable import parse_single_table # For VOTable parsing

import urllib.request # For URL handling
import urllib.parse # For URL parsing
from io import BytesIO # For byte handling
from glob import glob # For file handling

import matplotlib.pyplot as plt

# Colors for printing messages
RED = "\033[31m"       # For errors
GREEN = "\033[32;3m"   # For outputs
BLUE = "\033[34;1m"    # Request inputs
MAGENTA = "\033[35m"   # Header
CYAN = "\033[36m"      # Main words
ENDCOLOR = "\033[0m"

### RUN FOR ALL FILES IN DIRECTORY ###

INPUT_PATH = './data_output/'
OUTPUT_PATH = './check_figs/'

inputfiles_data = glob("%s/*_P3.fits" % INPUT_PATH)

# for i, inputfile_data in enumerate(inputfiles_data):

    # output_file = inputfile_data.replace(".fits", ".png").replace(INPUT_PATH, OUTPUT_PATH)
    # # output_file_data = inputfile_data.replace(".fits", "_data.png").replace(INPUT_PATH, OUTPUT_PATH)
    # # output_file_err = inputfile_data.replace(".fits", "_err.png").replace(INPUT_PATH, OUTPUT_PATH)

    # if inputfile_data == output_file:
    #     print(RED + "Error: input and output files are the same." + ENDCOLOR)
    #     sys.exit()

    # hdu_list = fits.open(inputfile_data)
    # hdu_data = hdu_list['DATA_EXT']
    # hdu_err = hdu_list['STAT_EXT']

    # sum_data = np.sum(hdu_data.data, axis=0)
    # max_data = np.max(hdu_data.data, axis=0)
    # sum_err = hdu_err.data[0]

    # fig = plt.figure(figsize=(10,10))
    # ax1 = fig.add_subplot(131)
    # ax2 = fig.add_subplot(132)
    # ax3 = fig.add_subplot(133)

    # ax1.imshow(sum_data, origin='lower')
    # ax2.imshow(max_data, origin='lower')
    # ax3.imshow(sum_err, origin='lower')

    # for ax in [ax1,ax2,ax3]:
    #     ax.set_yticklabels([])
    #     ax.set_xticklabels([])
    #     ax.grid(True, ls='--', color='w', alpha=0.3)

    # fig.tight_layout()
    # fig.savefig(output_file, bbox_inches='tight', dpi=300)

    # print('INITIALIZING PROCESS FOR FILE')
    # print("   Input file: ", inputfile_data)
    # print("   Output file: ", output_file)
    # print("\n")

# Define your percentile ranges for each row; adjust these as needed.
percentile_ranges = [(1, 99), (5, 95), (10, 90)]

for i, inputfile_data in enumerate(inputfiles_data):

    output_file = inputfile_data.replace(".fits", ".png").replace(INPUT_PATH, OUTPUT_PATH)

    if inputfile_data == output_file:
        print(RED + "Error: input and output files are the same." + ENDCOLOR)
        sys.exit()

    hdu_list = fits.open(inputfile_data)
    hdu_data = hdu_list['DATA_EXT']
    hdu_err = hdu_list['STAT_EXT']

    # Compute the images as before
    sum_data = np.sum(hdu_data.data, axis=0)
    max_data = np.max(hdu_data.data, axis=0)
    sum_err = hdu_err.data[0]

    # Prepare a 3x3 grid of subplots.
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    
    # Organize the three images and corresponding titles.
    data_images = [sum_data, max_data, sum_err]
    titles = ["Sum Data", "Max Data", "Sum Error"]

    # Loop over each row (each percentile range) and column (each image)
    for row, (p_low, p_high) in enumerate(percentile_ranges):
        for col, (image, title) in enumerate(zip(data_images, titles)):
            ax = axs[row, col]
            # Compute the vmin and vmax based on the given percentile range
            vmin = np.nanpercentile(image, p_low)
            vmax = np.nanpercentile(image, p_high)
            ax.imshow(image, origin='lower', vmin=vmin, vmax=vmax)
            ax.set_title(f"{title}\n{p_low}-{p_high} Percentile")
            # Remove tick labels
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            # Add a grid for visual reference
            ax.grid(True, ls='--', color='w', alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_file, bbox_inches='tight', dpi=300)

    print('INITIALIZING PROCESS FOR FILE')
    print("   Input file: ", inputfile_data)
    print("   Output file: ", output_file)
    print("\n")