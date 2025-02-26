"""
fits_duplicate_checker.py
-------------------------
This script processes FITS files in a specified directory to:
  - Pre-check each file to ensure the header 'PRODCATG' equals 'SCIENCE.CATALOGTILE'
  - Extract the 'ID' column from valid FITS files
  - Identify and report duplicate IDs across all processed files
  - Print informative messages for files that are skipped or contain duplicates

Example usage:
--------------
Run the script from the command line by specifying the input directory containing your FITS files. For example:

    python fits_duplicate_checker.py /diskb/phase3data/ftp/programs/KIDS/batch_29157/

This command will process all FITS files in the specified directory, pre-check them for the header condition, and then check for duplicate IDs.

Written by Ashley Barnes on 2025-02-26.
"""


import sys
import glob
import gc
from astropy.io import fits
import numpy as np

# Ensure an input directory was provided
if len(sys.argv) < 2:
    print("Usage: python script.py <input_directory>")
    sys.exit(1)

# Get the input directory from the command-line argument
input_dir = sys.argv[1]
if not input_dir.endswith('/'):
    input_dir += '/'

# Find all FITS files in the input directory
fits_files = glob.glob(input_dir + "*.fits")

print("\n--- FITS Duplicate Checker ---")

print("\nFound {} FITS files in directory {}".format(len(fits_files), input_dir))

# -------------------------------------------------------
# Pre-check files: only use files with header['PRODCATG'] == 'SCIENCE.CATALOGTILE'
# -------------------------------------------------------
valid_files = []
skipped_files = []
print("\nPre-checking FITS files for header 'PRODCATG' == 'SCIENCE.CATALOGTILE'.")
for filename in fits_files:
    try:
        with fits.open(filename, memmap=True) as hdul:
            header_value = hdul[0].header.get('PRODCATG')
            if header_value == 'SCIENCE.CATALOGTILE':
                valid_files.append(filename)
            else:
                skipped_files.append(filename)
                print("   Skipping file {}: header 'PRODCATG' value is '{}'".format(filename, header_value))
    except Exception as e:
        skipped_files.append(filename)
        print("   Skipping file {} due to error: {}".format(filename, e))

print("Total valid files after pre-check: {}".format(len(valid_files)))
if skipped_files:
    print("   Total skipped files: {}. Skipped files:".format(len(skipped_files)))
    for file in skipped_files:
        print("  {}".format(file))
else:
    print("   No files were skipped during pre-check.")

all_ids = []


# -------------------------------------------------------
# Process each valid FITS file and extract the "ID" column
# -------------------------------------------------------
print("\n---------------------------------")
print("Processing files")
for idx, filename in enumerate(valid_files, start=1):
    print("   File ({}/{}): {}".format(idx, len(valid_files), filename))
    try:
        with fits.open(filename, memmap=True) as hdul:
            try:
                table_data = hdul[1].data
            except IndexError:
                print("   No extension found in {}; skipping file.".format(filename))
                continue

            try:
                id_data = table_data['ID']
            except KeyError:
                print("   Column 'ID' not found in {}; skipping file.".format(filename))
                continue

            # Append the IDs array for later duplicate checks
            all_ids.append(id_data)
    except Exception as e:
        print("   Error processing {}: {}".format(filename, e))
    finally:
        # Try cleaning up to free memory
        try:
            del hdul, table_data, id_data
        except Exception:
            pass

    if idx % 100 == 0:
        gc.collect()
        print("[INFO] Performed garbage collection after processing {} files.".format(idx))

# -------------------------------------------------------
# Vectorized duplicate search: find duplicate IDs across all files
# -------------------------------------------------------
print("\n---------------------------------")
if all_ids:
    # Concatenate all ID arrays into one
    all_ids = np.concatenate(all_ids)
    print("Total number of IDs collected: {}".format(len(all_ids)))
    
    sorted_ids = np.sort(all_ids)
    dup_mask = sorted_ids[1:] == sorted_ids[:-1]
    duplicates_vector = sorted_ids[1:][dup_mask]
    duplicate_ids_vector = np.unique(duplicates_vector)
    
    if duplicate_ids_vector.size:
        print("\nDuplicate IDs found:")
        for dup in duplicate_ids_vector:
            print("   {}".format(dup))
else:
    print("No IDs were collected from the files.")

# -------------------------------------------------------
# Check each file for the presence of duplicate IDs
# -------------------------------------------------------
if duplicate_ids_vector.size:
    print("\nThe above duplicate IDs found in files:")
    for idx, filename in enumerate(valid_files, start=1):
        try:
            with fits.open(filename, memmap=True) as hdul:
                try:
                    table_data = hdul[1].data
                except IndexError:
                    print("   No extension found in {}; skipping file.".format(filename))
                    continue

                try:
                    id_data = table_data['ID']
                except KeyError:
                    print("   Column 'ID' not found in {}; skipping file.".format(filename))
                    continue

                for cid in duplicate_ids_vector:
                    if cid in id_data:
                        print("   {} found in file: {}".format(cid, filename))
        except Exception as e:
            print("   Error processing {}: {}".format(filename, e))
        finally:
            try:
                del hdul, table_data, id_data
            except Exception:
                pass

        if idx % 100 == 0:
            gc.collect()
else:
    print("\nNo duplicate IDs found in any file!")

print("\nScript execution complete -- Have a nice day!")