# **Spectrum** - Phase 3 Science Data Products (SDP) Submission Guide

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **spectral data**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

Spectral data products include one-dimensional spectra, extracted spectra from integral field units, and multi-object spectroscopy. These must follow strict formatting and metadata requirements to ensure archive compatibility.

## 2. Spectral Data Submission Requirements

### 2.1 General Data Format
- Spectral data must be provided in **FITS format** with appropriate metadata.
- The **spectrum data resides in the Primary HDU**.
- Supported formats:
  - **SCIENCE.SPECTRUM**: 1D extracted spectrum
  - **SCIENCE.CUBE.IFS**: Integral field unit spectral cube - see **CUBE.IFS** document 
- The **SPECSYS** keyword must specify the spectral reference system (e.g., 'BARYCENT', 'TOPOCENT').
- All spectra must include associated error estimates and metadata.

### 2.2 Wavelength Calibration
- Wavelength calibration must follow **FITS WCS conventions for spectral data**.
- **Wavelength coordinate system**:
  - Units: Angstrom (Å), nanometers (nm), or microns (µm)
  - Frequency (GHz) or energy (keV) are also acceptable for specific cases.
- **Spectral resolution and dispersion must be defined**:
  - `SPEC_BIN`: Spectral bin width.
  - `SPEC_RES`: Spectral resolution.
  - `WAVELMIN`, `WAVELMAX`: Minimum and maximum wavelengths covered.

### 2.3 Flux Calibration and Photometry Requirements
- **Flux calibration status** must be specified using `FLUXCAL`:
  - 'ABSOLUTE' if the spectrum is fully flux-calibrated.
  - 'UNCALIBRATED' if only relative flux calibration is applied.
- **Flux uncertainty quantification**:
  - `FLUXERR`: Error associated with each spectral bin.
  - `SPEC_ERR`: Spectral uncertainty.
- **Total flux and normalization**:
  - `TOT_FLUX`: Integrated flux value.
  - `CONTNORM`: Continuum normalization factor.

### 2.4 Temporal Information
- Observation time must be recorded using **Modified Julian Date (MJD)**.
- `MJD-END` must be greater than or equal to `MJD-OBS`.
- `EXPTIME` (exposure time per spectral bin) and `TEXPTIME` (total exposure time) must be included.

## 3. Required and Recommended FITS Header Keywords

#### 3.1.1 Primary HDU Keywords (Generic FITS Standard)

| Keyword  | Description | Type |
|----------|-------------|------|
| `SIMPLE`  | Conforms to FITS standard (must be `T`) | Mandatory |
| `BITPIX`  | Number of bits per pixel (e.g., `-32` for floating-point data) | Mandatory |
| `NAXIS`   | Number of data axes | Mandatory |
| `NAXIS1`  | Length of first data axis | Mandatory |
| `NAXIS2`  | Length of second data axis | Mandatory |
| `EXTEND`  | FITS file may contain extensions (`T` or `F`) | Mandatory |
| `DATE`    | Date of file creation | Recommended |
| `COMMENT` | Comment string | Optional |
| `HISTORY` | Processing history of the file | Optional |
| `CHECKSUM` | FITS checksum to verify data integrity | Mandatory |
| `DATASUM`  | Data checksum to verify HDU contents | Mandatory |

#### 3.1.2 Primary HDU Keywords

| Keyword  | Description | Type |
|----------|-------------|------|
| `PRODCATG` | Data product category | Mandatory |
| `ASSOCi` | Association identifier | Mandatory (31) |
| `ASSONi` | Associated file name | Mandatory (32) |
| `ASSOMi` | MD5 checksum of associated file | Mandatory (31) |
| `ORIGIN` | Institution responsible for data | Mandatory |
| `TELESCOP` | Telescope name | Mandatory |
| `INSTRUME` | Instrument name | Mandatory |
| `OBJECT` | Target name | Mandatory |
| `RA`, `DEC` | Right Ascension and Declination | Mandatory |
| `EQUINOX` | Equinox of coordinates | Mandatory (37) |
| `RADESYS` | Coordinate reference system | Mandatory |
| `TIMESYS` | Time system used | Mandatory (38) |
| `EXPTIME` | Exposure time per spectral bin | Mandatory |
| `TEXPTIME` | Total exposure time | Mandatory |
| `MJD-OBS` | Modified Julian Date of observation start | Mandatory |
| `MJD-END` | Modified Julian Date of observation end | Mandatory |
| `PROG_ID` | Program ID | Mandatory ESO |
| `OBIDi` | Observation ID | Mandatory ESO |
| `NCOMBINE` | Number of combined exposures | Mandatory ESO |
| `OBSTECH` | Observation technique | Mandatory |
| `FLUXCAL` | Flux calibration status | Mandatory |
| `PROCSOFT` | Data processing software | Mandatory |
| `REFERENC` | Bibliographic reference | Mandatory (40) |
| `PROVi` | Provenance information | Mandatory ESO |
| `BUNIT` | Physical units of flux | Not Allowed |
| `CDi_j` | WCS transformation matrix elements | Not Allowed |
| `SPECSYS` | Spectral reference system | Mandatory |
| `EXT_OBJ` | External object identifier | Mandatory (42) |
| `CONTNORM` | Continuum normalization factor | Mandatory |
| `TOT_FLUX` | Total flux value | Mandatory |
| `FLUXERR` | Flux uncertainty per bin | Mandatory (43) |
| `WAVELMIN` | Minimum wavelength coverage | Mandatory |
| `WAVELMAX` | Maximum wavelength coverage | Mandatory |
| `LAMRMS` | RMS deviation of wavelength solution | Optional |
| `LAMNLIN` | Number of fitted wavelength lines | Optional |
| `SPEC_BIN` | Spectral bin width | Mandatory |
| `SPEC_ERR` | Spectral error estimate | Recommended |
| `SPEC_SYE` | Systematic spectral error | Recommended |
| `RA_ERR` | RA uncertainty | Recommended |
| `DEC_ERR` | DEC uncertainty | Recommended |
| `STOKES` | Stokes parameter | Ask Support |
| `SNR` | Signal-to-noise ratio | Mandatory |
| `SPEC_RES` | Spectral resolution | Mandatory |
| `STREHL` | Strehl ratio (AO) | Mandatory (51) |
| `ARCFILE` | Archive file name | Reserved |
| `CHECKSUM` | FITS checksum | Mandatory |
| `DATASUM` | FITS data checksum | Mandatory |
| `ORIGFILE` | Original file name | Reserved |
| `P3ORIG` | Original Phase 3 identifier | Reserved |

#### **Notes:**

- **(31)** Mandatory when ancillary files are provided in a non-FITS format.
- **(32)** Mandatory when ancillary files are associated with the scientific data.
- **(37)** If RADESYS='FK5', then EQUINOX=2000.0 is mandatory.
- **(38)** Mandatory if the system used is other than UTC.
- **(40)** If a refereed publication is not available at the time of the data release, this value can be left empty.
- **(42)** Mandatory for externally submitted data products.
- **(43)** Applies to FLUXCAL='ABSOLUTE'. If FLUXCAL='UNCALIBRATED', set FLUXERR=-1.
- **(51)** Mandatory for adaptive optics (AO) observations.

#### 3.2.1 Extension HDU Keywords (Generic FITS Standard)

| Keyword  | Description | Type |
|----------|-------------|------|
| `XTENSION` | Type of FITS extension (`IMAGE`, `TABLE`, `BINTABLE`) | Mandatory |
| `BITPIX`   | Number of bits per pixel in the extension | Mandatory |
| `INHERIT`  | Inherited Primary HDU | Optional |
| `NAXIS`    | Number of data axes in the extension | Mandatory |
| `NAXIS1`   | Length of first data axis in the extension | Mandatory |
| `NAXIS2`   | Length of second data axis in the extension | Mandator |
| `PCOUNT`   | Parameter count (typically `0` for no varying arrays) | Mandatory (`BINTABLE` only) |
| `GCOUNT`   | Group count (typically `1` for standard FITS files) | Mandatory  (`BINTABLE` only)|
| `EXTNAME`  | Name of the extension HDU | Recommended |
| `CHECKSUM` | FITS checksum to verify data integrity | Mandatory |
| `DATASUM`  | Data checksum to verify extension contents | Mandatory |

#### 3.2.2 Extension HDU Keywords

| Keyword  | Description | Type |
|----------|-------------|------|
| `RA`, `DEC` | Right Ascension and Declination | Mandatory |
| `NELEM` | Number of elements in the data array | Mandatory |
| `VOCLASS` | Virtual Observatory classification | Mandatory |
| `VOPUB` | Virtual Observatory publication status | Mandatory |
| `TITLE` | Title of the spectrum | Mandatory |
| `APERTURE` | Aperture size used during observation | Mandatory |
| `TELAPSE` | Elapsed observation time | Mandatory |
| `TMID` | Midpoint time of observation | Mandatory |
| `SPEC_VAL` | Central spectral value | Mandatory |
| `SPEC_BW` | Spectral bandwidth | Mandatory |
| `ARCFILE` | Archive file name | Reserved |
| `CHECKSUM` | FITS checksum | Mandatory |
| `DATASUM` | FITS data checksum | Mandatory |
| `ORIGFILE` | Original file name | Reserved |
| `P3ORIG` | Original Phase 3 identifier | Reserved |

## Additional Considerations

- **Consistency Across HDUs:**  
  Keywords such as `RA`, `DEC`, `CHECKSUM`, and `DATASUM` must have identical values in both the primary and extension headers to maintain consistency.

- **MEF File Structure:**  
  For Multi-Extension FITS (MEF) files, include keywords like `SCIDATA`, `ERRDATA`, and `QUALDATA` when your dataset includes additional layers (e.g., separate error arrays or quality maps).

- **Validation:**  
  Prior to submission, run your FITS files through the ESO-provided validation tools to ensure that all header keywords comply with Phase 3 SDP standards.

For further details and examples of FITS headers, please refer to the [ESO Phase 3 FAQ](https://www.eso.org/sci/observing/phase3/faq.html) and the [ESO Phase 3 Overview](https://www.eso.org/sci/observing/phase3/overview.html).
## 4. File Structure and Validation
- **Single Spectrum FITS Format**: Data is stored in the primary HDU.
- **Multi-Extension FITS (MEF) Format**:
  - Science data in separate HDU.
  - Error estimates, variance, and quality maps stored in additional HDUs.
- Keywords for MEF files:
  - `SCIDATA`: Science data extension name
  - `ERRDATA`: Error data extension name
  - `QUALDATA`: Quality map extension name

## 5. Summary of Best Practices
- Follow the **FITS keyword requirements** for wavelength calibration, flux calibration, and data quality.
- Use **Multi-Extension FITS (MEF)** for spectral data with multiple layers.
- Ensure metadata consistency and validation before submission.
- Utilize **ESO-provided validation tools** to check compliance.

For further details and example FITS headers, refer to the [ESO Phase 3 FAQ](https://www.eso.org/sci/observing/phase3/faq.html) and the [ESO Phase 3 Overview](https://www.eso.org/sci/observing/phase3/overview.html).

