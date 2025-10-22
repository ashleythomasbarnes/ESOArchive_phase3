####
# Script to check the data cubes from the Phase 3 data.
# This script will create a 3x3 grid of images for each data cube, showing the sum of the data, the maximum of the data, and the sum of the error.
# The images will be saved in the same directory as the input data cube, with the same name but with a .png extension.
###

import sys
import numpy as np
from astropy.io import fits
from glob import glob
import matplotlib.pyplot as plt

# Colors for printing messages
RED = "\033[31m"       # For errors
GREEN = "\033[32;3m"   # For outputs
BLUE = "\033[34;1m"    # For request inputs
MAGENTA = "\033[35m"   # Header
CYAN = "\033[36m"      # Main words
ENDCOLOR = "\033[0m"

# INPUT PATH
INPUT_PATH = './checks/'
OUTPUT_PATH = './checks_figs/'
inputfiles_data = glob("%s/*.fits" % INPUT_PATH)

# Define the percentile ranges for the three columns.
percentile_ranges = [(1, 99), (5, 95), (15, 85)]

for i, inputfile_data in enumerate(inputfiles_data):
    print('PROCESSING FILE %d OF %d' % (i+1, len(inputfiles_data)))
    output_file = inputfile_data.replace(".fits", ".png").replace(INPUT_PATH, OUTPUT_PATH)
    
    if inputfile_data == output_file:
        print(RED + "Error: input and output files are the same." + ENDCOLOR)
        sys.exit()

    # Open the FITS file and extract the image data.
    hdu = fits.open(inputfile_data)[0]
    data = hdu.data

    # Create a 1x3 grid of subplots.
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    # Loop over the percentile ranges, applying each to the image.
    for j, (p_low, p_high) in enumerate(percentile_ranges):
        vmin = np.nanpercentile(data, p_low)
        vmax = np.nanpercentile(data, p_high)
        ax = axs[j]
        im = ax.imshow(data, origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"{p_low}-{p_high} Percentile")
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(True, ls='--', color='w', alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_file, bbox_inches='tight', dpi=300)
    print(GREEN + "Saved: " + output_file + ENDCOLOR)