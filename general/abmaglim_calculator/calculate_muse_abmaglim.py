#!/usr/bin/env python3
"""Calculate MUSE Phase 3 PIXNOISE and ABMAGLIM header values."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import convolve


# Constants and defaults used by the MUSE 2.8.x IDP and HDRL 1.4 code.
MUSE_FLUX_SCALE = 1e-20
VARIANCE_MIN = 0.0
VARIANCE_MAX = 25.0
VARIANCE_BIN_SIZE = 0.04
AB_ZERO_JY = 3631.0
FLAMBDA_TO_JY = 3.34e4
CPL_STD_MAD = 1.4826
CPL_FWHM_SIG = 2.354820045030949382
PIXNOISE_COMMENT = (
    "[erg.s**(-1).cm**(-2).angstrom**(-1)] pixel-to-pixel noise"
)
PIXNOISE_FITS_COMMENT = (
    "[erg.s**(-1).cm**(-2).angstrom**(-1)] pixel-to-"
)
ABMAGLIM_COMMENT = "5-sigma magnitude limit for point sources"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate pipeline-style PIXNOISE and ABMAGLIM values for a "
            "MUSE cube."
        )
    )
    parser.add_argument("cube", type=Path, help="MUSE FITS cube")
    return parser.parse_args()


def flux_scale_from_bunit(bunit: str | None) -> float:
    """Extract the effective power-of-ten scale from a MUSE BUNIT."""
    if not bunit:
        raise ValueError("The DATA extension is missing BUNIT.")

    # Variance units are commonly written either with an explicit 10**(-40)
    # scale or by squaring the complete 10**(-20) flux unit.
    compact_bunit = re.sub(r"\s+", "", bunit)
    squared_unit_match = re.match(
        r"^\(10\*\*\(([+-]?\d+)\).+\)\*\*\(?([+-]?\d+)\)?$",
        compact_bunit,
        flags=re.IGNORECASE,
    )
    if squared_unit_match:
        exponent, outer_power = map(int, squared_unit_match.groups())
        return 10.0 ** (exponent * outer_power)

    match = re.search(r"10\s*\*\*\s*\(\s*([+-]?\d+)\s*\)", bunit)
    if not match:
        raise ValueError(f"Cannot read the flux scale from BUNIT={bunit!r}.")

    return 10.0 ** int(match.group(1))


def spatial_pixel_scales_arcsec(header: fits.Header) -> tuple[float, float]:
    """Return the two spatial pixel scales in arcsec per pixel."""
    if "CD1_1" in header and "CD2_2" in header:
        xscale = math.hypot(header["CD1_1"], header.get("CD2_1", 0.0))
        yscale = math.hypot(header.get("CD1_2", 0.0), header["CD2_2"])
    elif "CDELT1" in header and "CDELT2" in header:
        xscale = abs(float(header["CDELT1"]))
        yscale = abs(float(header["CDELT2"]))
    else:
        raise ValueError(
            "The DATA extension needs a spatial CD matrix or CDELT1/CDELT2."
        )

    xscale *= 3600.0
    yscale *= 3600.0
    if xscale <= 0.0 or yscale <= 0.0:
        raise ValueError("The spatial pixel scales must be positive.")
    return xscale, yscale


def validate_inputs(
    cube_hdul: fits.HDUList,
) -> tuple[np.ndarray, np.ndarray, fits.Header, fits.Header]:
    """Validate the MUSE cube and return its relevant components."""
    try:
        cube_data_hdu = cube_hdul["DATA"]
        cube_stat_hdu = cube_hdul["STAT"]
    except KeyError as exc:
        raise ValueError("Expected DATA and STAT extensions in the cube.") from exc

    cube_data = cube_data_hdu.data
    cube_stat = cube_stat_hdu.data
    if cube_data is None or cube_stat is None:
        raise ValueError("One or more required cube extensions contain no data.")
    if cube_data.ndim != 3 or cube_stat.ndim != 3:
        raise ValueError("The cube DATA and STAT extensions must be three-dimensional.")
    if cube_data.shape != cube_stat.shape:
        raise ValueError("The cube DATA and STAT extension shapes do not match.")

    cube_scale = flux_scale_from_bunit(cube_data_hdu.header.get("BUNIT"))
    if cube_scale != MUSE_FLUX_SCALE:
        raise ValueError(
            "This MUSE implementation expects DATA values scaled by 10**(-20) "
            "erg s**(-1) cm**(-2) angstrom**(-1)."
        )

    stat_scale = flux_scale_from_bunit(cube_stat_hdu.header.get("BUNIT"))
    if stat_scale != MUSE_FLUX_SCALE**2:
        raise ValueError(
            "The STAT extension must contain variance scaled by 10**(-40)."
        )

    spatial_pixel_scales_arcsec(cube_data_hdu.header)

    return (
        cube_data,
        cube_stat,
        cube_hdul[0].header,
        cube_data_hdu.header,
    )


def compute_pixnoise(stat_cube: np.ndarray) -> float:
    """Reproduce the MUSE fixed-bin variance-mode PIXNOISE calculation."""
    nbins = int((VARIANCE_MAX - VARIANCE_MIN) / VARIANCE_BIN_SIZE) + 1
    histogram = np.zeros(nbins, dtype=np.int64)

    # Work one wavelength plane at a time so a multi-gigabyte cube remains
    # memory mapped rather than being copied into memory.
    for plane in stat_cube:
        valid = (
            np.isfinite(plane)
            & (plane >= VARIANCE_MIN)
            & (plane <= VARIANCE_MAX)
        )
        values = plane[valid]
        indices = (
            (values + VARIANCE_MIN) / VARIANCE_BIN_SIZE
        ).astype(np.int64)
        histogram += np.bincount(indices, minlength=nbins)

    peak_bin = int(np.argmax(histogram))
    mode_variance = (
        VARIANCE_MIN + (peak_bin + 0.5) * VARIANCE_BIN_SIZE
    )
    return math.sqrt(mode_variance) * MUSE_FLUX_SCALE


def collapse_cube_to_mean_image(data_cube: np.ndarray) -> np.ndarray:
    """Make the unfiltered finite-voxel mean image used by MUSE IDP."""
    total = np.zeros(data_cube.shape[1:], dtype=np.float64)
    count = np.zeros(data_cube.shape[1:], dtype=np.int32)

    # Plane-wise accumulation preserves the wavelength order used by the C
    # implementation and needs only a few spatial-size working arrays.
    for plane in data_cube:
        valid = np.isfinite(plane)
        np.add(total, plane, out=total, where=valid)
        count += valid

    collapsed = np.full(data_cube.shape[1:], np.nan, dtype=np.float32)
    np.divide(
        total,
        count,
        out=collapsed,
        where=count > 0,
        casting="unsafe",
    )
    return collapsed


def hdrl_gaussian_kernel(fwhm_pixels: float, size: int) -> np.ndarray:
    """Create the peak-normalized Gaussian kernel used by HDRL."""
    coordinate = np.arange(size, dtype=np.float64) - 0.5 * (size - 1)
    x, y = np.meshgrid(coordinate, coordinate)
    sigma = fwhm_pixels / math.sqrt(4.0 * math.log(4.0))
    return np.exp(-(x * x + y * y) / (2.0 * sigma * sigma))


def hdrl_convolve(image: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Apply HDRL-style reflected convolution while ignoring invalid pixels."""
    valid = np.isfinite(image)
    numerator = convolve(
        np.where(valid, image, 0.0),
        kernel,
        mode="reflect",
    )
    denominator = convolve(
        valid.astype(np.float64),
        kernel,
        mode="reflect",
    )

    result = np.full(image.shape, np.nan, dtype=np.float64)
    np.divide(numerator, denominator, out=result, where=denominator > 0.0)
    return result


def median_absolute_deviation(values: np.ndarray) -> float:
    median = float(np.median(values))
    return float(np.median(np.abs(values - median)))


def hdrl_weighted_mode(values: np.ndarray) -> float:
    """Reproduce HDRL's automatic-bin weighted histogram mode."""
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("No finite pixels are available for the HDRL mode.")

    mad = median_absolute_deviation(values)
    bin_size = (
        2.0
        * 3.49
        * (mad * CPL_STD_MAD)
        / values.size ** (1.0 / 3.0)
    )
    if bin_size <= 0.0:
        bin_size = np.nextafter(0.0, 1.0)

    histogram_min = float(np.min(values)) - 0.5 * bin_size
    initial_max = float(np.max(values)) + 0.5 * bin_size
    nbins = math.floor((initial_max - histogram_min) / bin_size) + 1
    histogram_max = histogram_min + nbins * bin_size

    histogram, _ = np.histogram(
        values,
        bins=nbins,
        range=(histogram_min, histogram_max),
    )
    peak_bin = int(np.argmax(histogram))
    peak_count = int(histogram[peak_bin])

    # HDRL averages the lower edges if several bins share the maximum.
    peak_bins = np.flatnonzero(histogram == peak_count)
    peak_level = float(np.mean(histogram_min + peak_bins * bin_size))

    previous_count = int(histogram[peak_bin - 1]) if peak_bin > 0 else 0
    next_count = (
        int(histogram[peak_bin + 1]) if peak_bin + 1 < nbins else 0
    )
    difference_left = peak_count - previous_count
    difference_right = peak_count - next_count
    denominator = difference_left + difference_right
    factor = difference_left / denominator if denominator else 0.5
    if factor == 0.0 or not math.isfinite(factor):
        factor = 0.5

    return peak_level + bin_size * factor


def compute_abmaglim(
    data_cube: np.ndarray,
    primary_header: fits.Header,
    data_header: fits.Header,
) -> float:
    """Reproduce the MUSE IDP call to HDRL's limiting-magnitude code."""
    for keyword in ("WAVELMIN", "WAVELMAX", "SKY_RES"):
        if keyword not in primary_header:
            raise ValueError(f"The cube primary header is missing {keyword}.")

    wavelength_min = float(primary_header["WAVELMIN"]) * 10.0
    wavelength_max = float(primary_header["WAVELMAX"]) * 10.0
    if wavelength_min <= 0.0 or wavelength_max <= wavelength_min:
        raise ValueError("WAVELMIN/WAVELMAX do not define a valid wavelength range.")

    xscale, yscale = spatial_pixel_scales_arcsec(data_header)
    fwhm_pixels = abs(float(primary_header["SKY_RES"])) / max(xscale, yscale)
    if fwhm_pixels <= 0.0:
        raise ValueError("SKY_RES must be non-zero.")

    kernel_size = int(fwhm_pixels) * 2 + 1
    collapsed = collapse_cube_to_mean_image(data_cube)
    kernel = hdrl_gaussian_kernel(fwhm_pixels, kernel_size)
    smoothed = hdrl_convolve(collapsed.astype(np.float64), kernel)
    finite_smoothed = smoothed[np.isfinite(smoothed)]

    mode = hdrl_weighted_mode(finite_smoothed)
    background = finite_smoothed[finite_smoothed <= mode]
    if background.size < 2:
        raise ValueError("Too few background pixels remain for ABMAGLIM.")

    mad = median_absolute_deviation(background)
    if mad <= 0.0:
        mad = np.nextafter(0.0, 1.0)

    correction = 1.0 / math.sqrt(1.0 - 2.0 / math.pi)
    noise = mad * CPL_STD_MAD * correction
    point_source_norm = (
        4.0 * math.pi * (fwhm_pixels / CPL_FWHM_SIG) ** 2
    )
    zeropoint = (
        -2.5
        * math.log10(
            FLAMBDA_TO_JY
            * wavelength_min
            * wavelength_max
            / AB_ZERO_JY
        )
        - 2.5 * math.log10(MUSE_FLUX_SCALE)
    )
    return (
        -2.5 * math.log10(5.0 * noise * point_source_norm)
        + zeropoint
    )


def calculate_muse_values(cube_path: Path) -> tuple[float, float]:
    """Calculate PIXNOISE and ABMAGLIM from a MUSE cube."""
    with fits.open(cube_path, memmap=True) as cube_hdul:
        cube_data, cube_stat, primary_header, data_header = validate_inputs(
            cube_hdul,
        )
        pixnoise = compute_pixnoise(cube_stat)
        abmaglim = compute_abmaglim(
            cube_data,
            primary_header,
            data_header,
        )
    return pixnoise, abmaglim


def format_header_values(pixnoise: float, abmaglim: float) -> tuple[str, str]:
    """Format the two values as Phase 3 header-style output lines."""
    return (
        f"PIXNOISE= {pixnoise:.14E} / {PIXNOISE_COMMENT}",
        f"ABMAGLIM= {abmaglim:.13f} / {ABMAGLIM_COMMENT}",
    )


def main() -> None:
    arguments = parse_arguments()
    if not arguments.cube.is_file():
        raise SystemExit(f"Cube not found: {arguments.cube}")

    try:
        pixnoise, abmaglim = calculate_muse_values(arguments.cube)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}") from exc

    for line in format_header_values(pixnoise, abmaglim):
        print(line)


if __name__ == "__main__":
    main()
