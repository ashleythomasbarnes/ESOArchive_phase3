import glob
from astropy.io import fits
import numpy as np
import gc
import os

# Define directories for each batch
# batch_29157_dir = "/diskb/phase3data/ftp/programs/KIDS/batch_29157/"
batch_29157_dir = "/diskb/phase3data/ftp/programs/KIDS/batch_28735/" # TEST
batch_28735_dir = "/diskb/phase3data/ftp/programs/KIDS/batch_28735/"

# Files in the new batch (batch_29157)
new_files = [
    "KiDS_DR5.0_1.2_-35.1_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_237.0_0.5_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_351.8_-31.2_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_222.6_2.5_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_28.0_-31.2_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_355.3_-32.1_ugri1i2ZYJHKs_cat.fits"
]

# Files in the old batch (batch_28735)
old_files = [
    "KiDS_DR5.0_2.4_-35.1_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_237.0_-0.5_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_221.6_2.5_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_26.8_-31.2_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_350.7_-31.2_ugri1i2ZYJHKs_cat.fits",
    "KiDS_DR5.0_355.3_-31.2_ugri1i2ZYJHKs_cat.fits"]

# Construct full paths for all files
new_file_paths = [os.path.join(batch_29157_dir, f) for f in new_files]
old_file_paths = [os.path.join(batch_28735_dir, f) for f in old_files]

# Combine all valid file paths into one list
valid_files = new_file_paths + old_file_paths
print("Processing {} files in total.".format(len(valid_files)))

all_ids = []

# -------------------------------------------------------
# Load files into memory and get IDs from each file
# -------------------------------------------------------
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
# Vectorized Sorting Approach - get duplicates across all files
# -------------------------------------------------------
if all_ids:
    # Concatenate all arrays of IDs into one array
    all_ids = np.concatenate(all_ids)
    print("\nTotal number of IDs collected: {}".format(len(all_ids)))
    
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

# -------------------------------------------------------
# Check for specific duplicate IDs in each file
# -------------------------------------------------------
print("\nDuplicate IDs found in each file:")
for idx, filename in enumerate(valid_files, start=1):
    try:
        with fits.open(filename, memmap=True) as hdul:
            try:
                table_data = hdul[1].data
            except IndexError:
                print("  No extension found in {}; skipping file.".format(filename))
                continue

            try:
                id_data = table_data['ID']
            except KeyError:
                print("  Column 'ID' not found in {}; skipping file.".format(filename))
                continue

            # Check if any of the duplicate IDs exist in this file
            for cid in duplicate_ids_vector:
                if cid in id_data:
                    print("  Duplicate ID {} found in file: {}".format(cid, filename))
                    
    except Exception as e:
        print("  Error processing {}: {}".format(filename, e))

    # Clean up references to free memory
    del hdul, table_data, id_data

    if idx % 100 == 0:
        gc.collect()
        #print("Performed garbage collection after processing {} files.".format(idx))

print("\nScript execution complete.")
