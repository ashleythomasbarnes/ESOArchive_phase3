#!/usr/bin/env python3
import sys
from glob import glob
import os
import tarfile
import hashlib
from astropy.io import fits


# ----------------------------
# Configuration
# ----------------------------
PRIMARY_SUFFIX = "_S1D_FINAL_A.fits"

ASSOCIATED_SUFFIXES = [
    "_S1D_FINAL_B.fits",
    "_CCF_A.fits",
    "_DRIFT_MATRIX_B.fits",
    "_S2D_A.fits",
    "_S2D_BLAZE_A.fits",
]

TAR_MEMBERS = ["_S2D_A.fits", "_S2D_BLAZE_A.fits"]
ASSOC5_VALUE = "ANCILLARY.2DECHELLE.TAR"


def md5sum(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def main(data_dir):
    print(f"Working directory: {data_dir}")
    assert os.path.isdir(data_dir), f"Not a directory: {data_dir}"

    # Always work with absolute paths internally
    data_dir = os.path.abspath(data_dir)

    primary_files = sorted(
        glob(os.path.join(data_dir, f"r.*{PRIMARY_SUFFIX}"))
    )
    print(f"Found {len(primary_files)} primary files")

    prefixes = [
        p[:-len(PRIMARY_SUFFIX)] for p in primary_files
    ]
    n = len(prefixes)

    # Count checks
    for suf in ASSOCIATED_SUFFIXES:
        files = glob(os.path.join(data_dir, f"r.*{suf}"))
        print(f"Found {len(files)} files matching *{suf}")
        assert len(files) == n, f"Count mismatch for {suf}"

    print("Count checks passed.\n")

    for prefix in prefixes:
        primary = prefix + PRIMARY_SUFFIX

        assoc_b = prefix + "_S1D_FINAL_B.fits"
        assoc_ccf = prefix + "_CCF_A.fits"
        assoc_drift = prefix + "_DRIFT_MATRIX_B.fits"

        s2d_a = prefix + "_S2D_A.fits"
        s2d_blaze = prefix + "_S2D_BLAZE_A.fits"

        tar_name = prefix + "_S1D_FINAL_A.tar"

        print(f"Processing {os.path.basename(prefix)}")

        # Existence checks
        for path in [primary, assoc_b, assoc_ccf, assoc_drift, s2d_a, s2d_blaze]:
            assert os.path.exists(path), f"Missing file: {path}"

        # Create tar
        print(f"  Creating {os.path.basename(tar_name)}")
        with tarfile.open(tar_name, "w") as tf:
            tf.add(s2d_a, arcname=os.path.basename(s2d_a))
            tf.add(s2d_blaze, arcname=os.path.basename(s2d_blaze))

        # Remove packed files
        print("  Removing S2D files")
        os.remove(s2d_a)
        os.remove(s2d_blaze)

        # MD5
        tar_md5 = md5sum(tar_name)

        # Header update
        print("  Updating FITS header")
        with fits.open(primary, mode="update") as hdul:
            hdr = hdul[0].header

            hdr["ASSON1"] = (os.path.basename(assoc_b), "Asson 1")
            hdr["ASSON2"] = (os.path.basename(assoc_ccf), "Asson 2")
            hdr["ASSON3"] = (os.path.basename(assoc_drift), "Asson 3")
            hdr["ASSON4"] = (os.path.basename(tar_name), "Asson tar")

            hdr["ASSOC4"] = (ASSOC5_VALUE, "Asson class")
            hdr["ASSOM4"] = (tar_md5, "MD5 checksum")

            # Update FITS integrity keywords for the primary HDU
            hdul[0].add_datasum()
            hdul[0].add_checksum()

            hdul.flush()

        print("  Done.\n")

    print("All files processed successfully.")


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Usage: python3 update_headers.py <data_dir>"
    main(sys.argv[1])