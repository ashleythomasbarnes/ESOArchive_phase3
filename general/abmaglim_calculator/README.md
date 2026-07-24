# MUSE ABMAGLIM and PIXNOISE calculator

This tool calculates two ESO Phase 3 primary-header values directly from a
MUSE FITS cube:

- `PIXNOISE`: the typical noise in one spatial and wavelength pixel.
- `ABMAGLIM`: the faintest point source expected to be detected at a
  signal-to-noise ratio of 5.

The calculation follows the MUSE 2.8.x IDP code and the HDRL 1.4
implementation.

## Requirements

The input must be a three-dimensional MUSE cube containing:

- a `DATA` extension with flux values in units scaled by `10**(-20)`;
- a matching `STAT` extension with variance values in units scaled by
  `10**(-40)`;
- the primary-header keywords `WAVELMIN`, `WAVELMAX`, and `SKY_RES`; and
- spatial pixel-scale information in the `DATA` header, supplied through a CD
  matrix or `CDELT1` and `CDELT2`.

The script checks these items before starting the calculation and reports an
error if a required extension, keyword, unit, or image dimension is missing or
inconsistent.

## Install

Install the three Python dependencies, or use an environment that already has
Astropy, NumPy, and SciPy:

```bash
python -m pip install -r requirements.txt
```

## Calculate the values

From this directory, run:

```bash
python calculate_muse_abmaglim.py \
  '/path/to/muse_cube.fits'
```

This command only reads the cube. It does not change the file. The output has
the following form:

```text
PIXNOISE= 1.80554700852678E-20 / [erg.s**(-1).cm**(-2).angstrom**(-1)] pixel-to-pixel noise
ABMAGLIM= 24.4487056104046 / 5-sigma magnitude limit for point sources
```

The numbers above are an example. The values for another cube will normally be
different.

To calculate the values and add or replace the `PIXNOISE` and `ABMAGLIM`
keywords in the cube's primary header, run:

```bash
python calculate_replace_muse_header_values.py \
  '/path/to/writable_muse_cube.fits'
```

This second command changes the input FITS file in place and does not create a
backup. Make a copy first if the original file must be preserved. The script
prints the values only after reopening the primary header and confirming that
both values were stored. If the primary HDU already contains FITS checksum
cards, the script refreshes its checksum.

## What the calculation does

The two values are calculated separately. `PIXNOISE` comes from the `STAT`
variance cube. `ABMAGLIM` comes from the science values in the `DATA` cube and
does not use the calculated `PIXNOISE` value.

### 1. Check the cube

The script first confirms that `DATA` and `STAT` are both three-dimensional
and have the same shape. It also checks their units, the wavelength limits, the
spatial resolution, and the spatial pixel scale. Invalid values such as `NaN`
and infinity are ignored during the numerical calculations.

### 2. Calculate `PIXNOISE`

The `STAT` extension stores a variance for every voxel. A voxel is one element
of the cube, at one sky position and one wavelength. Variance is noise squared,
so the script must estimate a representative variance and then take its square
root.

To match the MUSE pipeline, the script:

1. Keeps finite `STAT` values between 0 and 25.
2. Sorts those values into fixed-width histogram bins of 0.04.
3. Finds the bin containing the largest number of values. This represents the
   most common, or typical, variance in the cube.
4. Uses the center of that bin as the representative variance.
5. Takes its square root and applies the MUSE `10**(-20)` flux scale.

The result is `PIXNOISE`, expressed as a typical pixel-to-pixel flux-density
noise. This is a pipeline-wide estimate from the variance cube, not a
measurement from a hand-selected blank-sky region.

### 3. Collapse the `DATA` cube into one image

For every sky pixel, the script averages all finite values along the wavelength
axis. This produces a two-dimensional, unfiltered mean image while preserving
the cube's flux-density units.

No external white-light image is used. An archived `white` filter image is not
numerically identical to this internal unfiltered mean image, so substituting
one would not reproduce the pipeline header value.

### 4. Convert the image resolution into pixels

`SKY_RES` gives the spatial resolution, or seeing full width at half maximum
(FWHM), in arcseconds. The script reads the angular size of a spatial pixel
from the `DATA` header and divides `SKY_RES` by the larger of the two spatial
pixel scales. The result is the FWHM in pixels.

In simple terms, this tells the calculation how many image pixels a point
source covers.

### 5. Smooth the image as a point source would appear

The mean image is smoothed with a Gaussian kernel whose FWHM matches the value
calculated in the previous step. This combines neighboring pixels over the
area occupied by a point source. At image edges and around invalid pixels, the
calculation adjusts the normalization so that missing values do not
artificially lower the result.

### 6. Estimate the background noise

The script makes an automatically sized histogram of all finite pixels in the
smoothed image and estimates its mode, meaning the most common background
level. It then keeps only pixels at or below that level. This lower half is
used to reduce contamination from stars, galaxies, and other positive sources.

The spread of these background pixels is measured with the median absolute
deviation (MAD), a statistic that is less sensitive to unusual pixels than a
standard deviation. The pipeline correction factors convert this one-sided MAD
into an estimate of the full background noise.

### 7. Calculate a five-sigma point-source flux limit

The background noise is scaled by the effective area of a Gaussian point
source, using the FWHM in pixels. The result is then multiplied by 5. This is
the flux-density level at which a point source would have a signal-to-noise
ratio of 5 under the pipeline assumptions.

### 8. Convert the flux limit to an AB magnitude

The script uses `WAVELMIN` and `WAVELMAX` from the primary header to construct
the wavelength-dependent AB conversion. MUSE Phase 3 wavelength limits are
stored in nanometres, so the script converts them to angstroms. It also applies
the cube's `10**(-20)` flux scale and the AB reference flux of 3631 Jy.

Finally, the five-sigma flux-density limit is converted to `ABMAGLIM`. Because
the magnitude scale runs backwards, a larger `ABMAGLIM` means a fainter source
can be detected and therefore indicates a deeper observation.

## Practical interpretation

- `PIXNOISE` describes the typical noise of an individual cube voxel.
  Smaller values indicate less noise.
- `ABMAGLIM` describes the pipeline-estimated five-sigma detection limit for a
  point source. Larger magnitudes indicate greater depth.
- `ABMAGLIM` assumes a Gaussian point source and a representative background.
  It is not a completeness limit and does not guarantee that every source at
  that magnitude will be detected.
- These are pipeline-style global values for the cube. Noise and sensitivity
  can still vary across the field of view and with wavelength.

The cube is processed one wavelength plane at a time. The full 3D arrays are
kept memory-mapped, which avoids loading a multi-gigabyte MUSE cube into RAM.
