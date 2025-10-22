####
# Script to check the data cubes from the Phase 3 data.
# This script will create a 3x3 grid of images for each data cube, showing the sum of the data, the maximum of the data, and the sum of the error.
# The images will be saved in the same directory as the input data cube, with the same name but with a .png extension.
###

# Imports
import sys # For system handling
import numpy as np # For numerical operations
from astropy.io import fits # For FITS file handling
from glob import glob # For file handling
import matplotlib.pyplot as plt # For plotting

# Colors for printing messages
RED = "\033[31m"       # For errors
GREEN = "\033[32;3m"   # For outputs
BLUE = "\033[34;1m"    # Request inputs
MAGENTA = "\033[35m"   # Header
CYAN = "\033[36m"      # Main words
ENDCOLOR = "\033[0m"

### RUN FOR ALL FILES IN DIRECTORY ###

# INPUT
# INPUT_PATH = './data_output/'
INPUT_PATH = '/diskb/phase3data/ftp/programs/SEDIGISM/batch_29448/'

# OUTPUT 
OUTPUT_PATH = './checks/'

inputfiles_data = glob("%s/*_P3.fits" % INPUT_PATH)

# Define your percentile ranges for each row; adjust these as needed.
percentile_ranges = [(1, 99), (5, 95), (10, 90)]

# for i, inputfile_data in enumerate(inputfiles_data):

#     print('PROCESSING FILE %d OF %d' % (i+1, len(inputfiles_data)))

#     output_file = inputfile_data.replace(".fits", ".png").replace(INPUT_PATH, OUTPUT_PATH)

#     if inputfile_data == output_file:
#         print(RED + "Error: input and output files are the same." + ENDCOLOR)
#         sys.exit()

#     hdu_list = fits.open(inputfile_data)
#     hdu_data = hdu_list['DATA_EXT']
#     hdu_err = hdu_list['STAT_EXT']

#     # Compute the images as before
#     sum_data = np.sum(hdu_data.data, axis=0)
#     max_data = np.max(hdu_data.data, axis=0)
#     sum_err = hdu_err.data[0]

#     # Prepare a 3x3 grid of subplots.
#     fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    
#     # Organize the three images and corresponding titles.
#     data_images = [sum_data, max_data, sum_err]
#     titles = ["Sum Data", "Max Data", "Sum Error"]

#     # Loop over each row (each percentile range) and column (each image)
#     for row, (p_low, p_high) in enumerate(percentile_ranges):
#         for col, (image, title) in enumerate(zip(data_images, titles)):
#             ax = axs[row, col]
#             # Compute the vmin and vmax based on the given percentile range
#             vmin = np.nanpercentile(image, p_low)
#             vmax = np.nanpercentile(image, p_high)
#             ax.imshow(image, origin='lower', vmin=vmin, vmax=vmax)
#             ax.set_title(f"{title}\n{p_low}-{p_high} Percentile")
#             # Remove tick labels
#             ax.set_xticklabels([])
#             ax.set_yticklabels([])
#             # Add a grid for visual reference
#             ax.grid(True, ls='--', color='w', alpha=0.3)

#     fig.tight_layout()
#     fig.savefig(output_file, bbox_inches='tight', dpi=300)

#     print("   Input file: ", inputfile_data)
#     print("   Output file: ", output_file)
#     print("\n")

for i, inputfile_data in enumerate(inputfiles_data):

    print('PROCESSING FILE %d OF %d' % (i+1, len(inputfiles_data)))

    output_file_sum = inputfile_data.replace(".fits", "_sum.fits").replace(INPUT_PATH, OUTPUT_PATH)
    output_file_max = inputfile_data.replace(".fits", "_max.fits").replace(INPUT_PATH, OUTPUT_PATH)
    output_file_err = inputfile_data.replace(".fits", "_sum_err.fits").replace(INPUT_PATH, OUTPUT_PATH)

    if inputfile_data == output_file_sum or inputfile_data == output_file_max or inputfile_data == output_file_err:
        print(RED + "Error: input and output files are the same." + ENDCOLOR)
        sys.exit()

    hdu_list = fits.open(inputfile_data)
    # hdu_data = hdu_list['DATA_EXT']
    hdu_err = hdu_list['STAT_EXT']

    # # Compute the images as before
    # sum_data = np.sum(hdu_data.data, axis=0)
    # max_data = np.max(hdu_data.data, axis=0)
    sum_err = hdu_err.data[0]

    # hdu_data_sum = fits.PrimaryHDU(sum_data, header=hdu_data.header)
    # hdu_data_max = fits.PrimaryHDU(max_data, header=hdu_data.header)
    hdu_err_sum = fits.PrimaryHDU(sum_err, header=hdu_err.header)

    hdu_err_sum.writeto(output_file_err, overwrite=True)