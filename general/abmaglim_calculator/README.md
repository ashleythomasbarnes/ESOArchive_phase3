# MUSE ABMAGLIM and PIXNOISE calculator

`calculate_muse_abmaglim.py` calculates the two Phase 3 primary-header values
directly from a MUSE cube. The calculation follows the MUSE 2.8.x IDP code and
the HDRL 1.4 implementation supplied in `scripts_pipeline/`.

## Run

Install the three Python dependencies, or use an environment that already has
Astropy, NumPy, and SciPy:

```bash
python -m pip install -r requirements.txt
```

To calculate and print the values without modifying the cube:

```bash
python calculate_muse_abmaglim.py \
  'test_data/ADP.2025-08-01T13:23:18.841.fits'
```

Expected output for the included test data:

```text
PIXNOISE= 1.80554700852678E-20 / [erg.s**(-1).cm**(-2).angstrom**(-1)] pixel-to-pixel noise
ABMAGLIM= 24.4487056104046 / 5-sigma magnitude limit for point sources
```

To calculate the values and add or replace them in the cube's primary header:

```bash
python calculate_replace_muse_header_values.py \
  '/path/to/writable_muse_cube.fits'
```

The replacement script updates the input cube in place without creating a
backup. It prints the two values only after rereading the primary header and
verifying that both were stored successfully. If the primary HDU already has
FITS checksum cards, its checksum is refreshed after the update.

## What the calculation does

- `PIXNOISE` is the square root of the mode of the cube's `STAT` variance
  extension. MUSE uses a fixed histogram from 0 to 25 with bin size 0.04, then
  uses the center of the most populated bin.
- `ABMAGLIM` uses the MUSE IDP call to HDRL: a finite-voxel mean collapse of
  the cube, Gaussian smoothing based on `SKY_RES`, HDRL's weighted histogram
  mode, a one-sided background MAD estimate, and the five-sigma point-source
  normalization.
- No external whitelight image is needed. The archived `white` filter image is
  not numerically identical to the unfiltered mean image generated internally
  from the cube, so using that image directly would not reproduce the pipeline
  header value.

The cube is processed one wavelength plane at a time. The full 3D arrays are
kept memory-mapped, which avoids loading a multi-gigabyte MUSE cube into RAM.
