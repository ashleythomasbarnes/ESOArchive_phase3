import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
# from tqdm import tqdm
import glob
import warnings
warnings.filterwarnings("ignore") # Be CAREFUL...

# Function to list files matching a specified pattern
def listfiles(fname, path, keep_path=False):
    # Get all file paths matching the pattern using glob
    filenames = glob.glob('%s%s' % (path, fname))

    return filenames

# Define the base path where the search will be performed
inputpath = '/diskb/phase3data/ftp/programs/SEDIGISM/batch_22605/'

# Use the listfiles function to find all FITS files in the './data/' directory 
# that match the pattern '*_DR1.fits'. The keep_path flag ensures the full path is kept.
fnames = listfiles('*_DR1.fits', inputpath, keep_path=True)

# Set the number of line-free channels to use for noise estimation
nchans = 100

# Set S/N treshold cut
threshold = 2

# Loop through each FITS file with a progress bar provided by tqdm
for i, fname in enumerate(fnames):

    ####
    ## LETS TAKE THE BREAKS OFF THIS THING!
    # if i > 0: 
    #     continue
    ####

    print("Running file %i out of %i" %(i, len(fnames)))
    print("   "+fname)

    # Open the FITS file and get the primary HDU (Header/Data Unit)
    hdu = fits.open(fname)[0]

    # Remove any singleton dimensions from the data array (e.g., shape (1, n, m) becomes (n, m))
    hdu.data = hdu.data.squeeze()
    # Get the dimensions of the data cube (assumed to be 3D)
    sizes = hdu.data.shape
    
    # Select channels from the end of the cube (excluding the very last channel)
    free1 = hdu.data[sizes[0]-(nchans+1):sizes[0]-1,:,:]
    # Select channels from the beginning of the cube
    free2 = hdu.data[0:nchans,:,:]

    # Combine the two selected sets of channels along the channel axis
    free = np.concatenate([free2, free1], axis=0)

    # Calculate the noise (rms) using the median absolute deviation (MAD) method,
    # which is robust against outliers. The calculation is done across the channels.
    rmsmap = mad_std(free, axis=0, ignore_nan=True)

    # Create a boolean mask where the data exceeds twice the computed rms value
    mask = free > rmsmap * threshold
    # Make a copy of the "free" data to apply masking
    free_masked = free.copy()
    # Replace the high-value pixels (likely containing signals) with NaN
    free_masked[mask] = np.nan

    # Recalculate the rms noise after masking out high values to get a cleaner noise estimate
    rmsmap = mad_std(free_masked, axis=0, ignore_nan=True)

    # Expand the 2D rms map into a 3D cube by repeating it for each channel of the original data cube
    rmscube = np.repeat(rmsmap[np.newaxis, :, :], hdu.data.shape[0], axis=0)

    # Create a new Primary HDU with the rms cube and copy the original header
    rmscube = fits.PrimaryHDU(rmscube, hdu.header)
    # Write the rms cube to a new FITS file; the new file name appends '_RMSCUBE' to the original name.
    rmscube.writeto(fname.split('.fits')[0] + '_RMSCUBE.fits', overwrite=True)