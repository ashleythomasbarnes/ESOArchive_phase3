import glob
from astropy.io import fits
import numpy as np

# Input directory containing the FITS files
input_dir = "/diskb/phase3data/ftp/programs/KIDS/batch_28735/"

# Define the pattern for the filenames
file_pattern = input_dir + "*_ugri1i2ZYJHKs_cat.fits"

# Get a list of all FITS files matching the pattern
fits_files = glob.glob(file_pattern)
print("Found {} files matching the pattern.".format(len(fits_files)))

# Define the file to be excluded
excluded_file = input_dir + "KiDS_DR5.0_ugri1i2ZYJHKs_cat.fits"

# Filter out the excluded file
valid_files = [filename for filename in fits_files if filename != excluded_file]
if len(valid_files) < len(fits_files):
    print("Excluding file: {}".format(excluded_file))
print("Processing {} files after exclusion.".format(len(valid_files)))

all_ids = []

# Process each valid FITS file
for idx, filename in enumerate(valid_files, start=1):
    print("Processing file ({}/{}): {}".format(idx, len(valid_files), filename))
    try:
        # Open the FITS file with memory mapping enabled
        with fits.open(filename, memmap=True) as hdul:
            # Assume that the table is in the first extension (index 1)
            try:
                table_data = hdul[1].data
            except IndexError:
                print("  No extension found in {}; skipping file.".format(filename))
                continue

            # Extract the "ID" column; skip the file if the column is missing
            try:
                id_data = table_data['ID']
            except KeyError:
                print("  Column 'ID' not found in {}; skipping file.".format(filename))
                continue

            # Append the array of IDs to the list
            all_ids.append(id_data)
    except Exception as e:
        print("  Error processing {}: {}".format(filename, e))

# If we collected any IDs, concatenate them into one array and check for duplicates
if all_ids:
    # Concatenate all arrays of IDs into one array
    all_ids = np.concatenate(all_ids)
    print("Total number of IDs collected: {}".format(len(all_ids)))

    # Get the unique IDs using np.unique (without return_counts)
    unique_ids = np.unique(all_ids)
    # Manually count the number of occurrences for each unique ID
    counts = [np.sum(all_ids == uid) for uid in unique_ids]
    # Convert counts to a NumPy array so we can index with a boolean mask
    counts = np.array(counts)
    duplicate_ids = unique_ids[counts > 1]

    if duplicate_ids.size:
        print("Duplicate IDs found:")
        for dup in duplicate_ids:
            print(dup)
    else:
        print("No duplicate IDs found.")
else:
    print("No IDs were collected from the files.")
