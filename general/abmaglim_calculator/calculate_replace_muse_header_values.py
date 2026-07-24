#!/usr/bin/env python3
"""Calculate and replace MUSE Phase 3 PIXNOISE and ABMAGLIM values."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from astropy.io import fits

from calculate_muse_abmaglim import (
    ABMAGLIM_COMMENT,
    PIXNOISE_FITS_COMMENT,
    calculate_muse_values,
    format_header_values,
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate pipeline-style PIXNOISE and ABMAGLIM values and "
            "replace them in a MUSE cube's primary header."
        )
    )
    parser.add_argument("cube", type=Path, help="writable MUSE FITS cube")
    return parser.parse_args()


def write_header_values(
    cube_path: Path,
    pixnoise: float,
    abmaglim: float,
) -> None:
    """Add or replace the two values in the cube's primary header."""
    stored_pixnoise = float(f"{pixnoise:.14E}")
    stored_abmaglim = float(f"{abmaglim:.13f}")

    with fits.open(cube_path, mode="update", memmap=True) as cube_hdul:
        primary_hdu = cube_hdul[0]
        header = primary_hdu.header
        had_checksum = "CHECKSUM" in header or "DATASUM" in header

        header["PIXNOISE"] = (stored_pixnoise, PIXNOISE_FITS_COMMENT)
        header["ABMAGLIM"] = (stored_abmaglim, ABMAGLIM_COMMENT)

        if had_checksum:
            primary_hdu.add_checksum()

        cube_hdul.flush(output_verify="fix", verbose=False)


def verify_header_values(
    cube_path: Path,
    pixnoise: float,
    abmaglim: float,
) -> None:
    """Verify that the recalculated values were stored successfully."""
    header = fits.getheader(cube_path, 0)
    expected_pixnoise = float(f"{pixnoise:.14E}")
    expected_abmaglim = float(f"{abmaglim:.13f}")
    stored_pixnoise = header.get("PIXNOISE")
    stored_abmaglim = header.get("ABMAGLIM")
    if stored_pixnoise != expected_pixnoise:
        raise OSError("PIXNOISE verification failed after updating the cube.")
    if stored_abmaglim != expected_abmaglim:
        raise OSError("ABMAGLIM verification failed after updating the cube.")


def main() -> None:
    arguments = parse_arguments()
    if not arguments.cube.is_file():
        raise SystemExit(f"Cube not found: {arguments.cube}")
    if not os.access(arguments.cube, os.W_OK):
        raise SystemExit(f"Cube is not writable: {arguments.cube}")

    try:
        pixnoise, abmaglim = calculate_muse_values(arguments.cube)
        write_header_values(arguments.cube, pixnoise, abmaglim)
        verify_header_values(arguments.cube, pixnoise, abmaglim)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}") from exc

    for line in format_header_values(pixnoise, abmaglim):
        print(line)


if __name__ == "__main__":
    main()
