import glob
from astropy.io import fits
import numpy as np
from collections import Counter
import gc

# Input directory containing the FITS files
# input_dir = "/diskb/phase3data/ftp/programs/KIDS/batch_28735/"
input_dir = "/diskb/phase3data/ftp/programs/KIDS/batch_29157/"

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


# -------------------------------------------------------
# Load files into memory and get IDs from each file
# -------------------------------------------------------
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

    # Clean up references to free memory
    del hdul, table_data, id_data

    # Perform garbage collection every 100 files
    if idx % 100 == 0:
        gc.collect()
        print("Performed garbage collection after processing {} files.".format(idx))


# -------------------------------------------------------
# Vectorized Sorting Approach - get duplicates
# -------------------------------------------------------
# If we collected any IDs, concatenate them into one array and check for duplicates
if all_ids:
    # Concatenate all arrays of IDs into one array
    all_ids = np.concatenate(all_ids)
    print("Total number of IDs collected: {}".format(len(all_ids)))
    
    sorted_ids = np.sort(all_ids)
    # Compare adjacent elements (all but the first compared with all but the last)
    dup_mask = sorted_ids[1:] == sorted_ids[:-1]
    # Extract duplicates; they appear consecutively in the sorted array
    duplicates_vector = sorted_ids[1:][dup_mask]
    # Get unique duplicate IDs
    duplicate_ids_vector = np.unique(duplicates_vector)
    
    if duplicate_ids_vector.size:
        print("\nDuplicate IDs found (Vectorized Sorting Approach):")
        for dup in duplicate_ids_vector:
            print(dup)
    else:
        print("\nNo duplicate IDs found (Vectorized Sorting Approach).")



## -------------------------------------------------------
## Check for specific duplicate IDs in the files
## -------------------------------------------------------
print("\nDuplicate IDs found in files:")
# Process each valid FITS file
for idx, filename in enumerate(valid_files, start=1):
    # print("Processing file ({}/{}): {}".format(idx, len(valid_files), filename))
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

            # Check if any of the IDs in check_ids exist in this file
            for cid in duplicate_ids_vector:
                if cid in id_data:
                    print("  {} in file: {}".format(cid, filename))
                    
    except Exception as e:
        print("  Error processing {}: {}".format(filename, e))

    # Clean up references to free memory
    del hdul, table_data, id_data

    # Perform garbage collection every 100 files
    if idx % 100 == 0:
        gc.collect()
        # print("Performed garbage collection after processing {} files.".format(idx))

print("\nScript execution complete.")
