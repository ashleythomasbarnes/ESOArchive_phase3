# **Flux Maps** - Phase 3 Science Data Products (SDP) Submission Guide

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **flux maps**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

Flux maps are specialized imaging products where pixel values represent flux measurements in physical units (e.g. milli-Jansky per beam) - as, for example, observations using the bolometer arrays LABOCA and ARTEMIS at the APEX 12-metre telescope. They must follow specific formatting and metadata requirements to be valid for ESO archive ingestion.

## 2. Flux Map Submission Requirements

### 2.1 General Data Format
- Flux maps must be provided in **FITS format** with appropriate metadata.
- The **flux map data resides in the Primary HDU**.
- Supported formats:
  - **Single images**: Simple FITS (SCIENCE.IMAGE.FLUXMAP)
  - **Multi-extension images**: Multi-Extension FITS (MEF) with metadata extensions.
- The **BUNIT** keyword must specify the physical units of flux (e.g., Jy/pixel, W/m²/sr, or MJy/sr).
- All flux maps must include associated error and quality maps.

### 2.2 Astrometry Requirements
- World Coordinate System (WCS) must follow the **FITS convention**.
- **Celestial reference system**: International Celestial Reference System (ICRS).
- **Astrometric error quantification**:
  - `CRDER1`, `CRDER2`: Random error in each axis.
  - `CSYER1`, `CSYER2`: Systematic error in each axis.
- `RADESYS` must be set to either 'ICRS' or 'FK5'. If 'FK5' is used, `EQUINOX` must be `2000.0`.
- RA/DEC values must be in the correct range: `[0,360]`, `[-90,90]`.

### 2.3 Flux Calibration and Photometry Requirements
- **Flux calibration status** must be specified using `FLUXCAL`:
  - 'ABSOLUTE' if the flux calibration is performed.
  - 'UNCALIBRATED' if no absolute flux calibration is applied.
- The **physical flux units** should be set via `BUNIT`.
- **Flux uncertainty quantification**:
  - `FLUXERR`: Error associated with each pixel’s flux value.
- The **photometric zero point** must be defined with `PHOTZP`, if applicable.

### 2.4 Temporal Information
- Observation time must be recorded using **Modified Julian Date (MJD)**.
- `MJD-END` must be greater than or equal to `MJD-OBS`.
- `EXPTIME` (exposure time per pixel) and `TEXPTIME` (total exposure time) must be included.

## 3. Required and Recommended FITS Header Keywords

### 3.1 Primary HDU Keywords

Assuming the **flux map data resides in the Primary HDU** - see **Multi-extension images** for additional options. 

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
| `ASSONi` | Associated file name | Mandatory (33) |
| `ASSOMi` | MD5 checksum of associated file | Mandatory (31) |
| `ORIGIN` | Institution responsible for data (e.g. **APEX**) | Mandatory |
| `TELESCOP` | Telescope name | Mandatory |
| `INSTRUME` | Instrument name | Mandatory |
| `FILTER` | Filter used during observation | Mandatory |
| `OBJECT` | Target name | Mandatory |
| `RA`, `DEC` | Right Ascension and Declination | Mandatory |
| `EQUINOX` | Equinox of coordinates | Mandatory (37) |
| `RADESYS` | Coordinate reference system | Mandatory |
| `TIMESYS` | Time system used | Mandatory (38) |
| `EXPTIME` | Exposure time per pixel | Mandatory |
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
| `BUNIT` | Physical units of flux | Mandatory |
| `CRVALi` | Coordinate value at reference pixel | Mandatory |
| `CRPIXi` | Reference pixel coordinates | Mandatory |
| `CTYPEi` | Coordinate system type | Mandatory |
| `CUNITi` | Coordinate units | Mandatory |
| `CDi_j` | WCS transformation matrix elements | Mandatory |
| `CSYERi` | Systematic error | Recommended |
| `CRDERi` | Random error | Recommended |
| `FLUXERR` | Flux uncertainty per pixel | Mandatory |
| `WAVELMIN` | Minimum wavelength coverage | Mandatory |
| `WAVELMAX` | Maximum wavelength coverage | Mandatory |
| `RA_ERR` | RA uncertainty | Not Allowed |
| `DEC_ERR` | DEC uncertainty | Not Allowed |
| `BNOISE` | Background noise level | Mandatory |
| `MAPMODE` | Mapping mode used for observation | Mandatory |
| `FEBEi` | Backend system identifier | Mandatory |
| `STOKES` | Stokes parameter | Ask Support |
| `ABMAGLIM` | Absolute magnitude limit | Not Allowed |
| `SKY_RES` | Sky resolution | Mandatory |
| `SKY_RERR` | Sky resolution error | Mandatory (50) |
| `ARCFILE` | Archive file name | Reserved |
| `CHECKSUM` | FITS checksum | Mandatory |
| `DATASUM` | FITS data checksum | Mandatory |
| `ORIGFILE` | Original file name | Reserved |
| `P3ORIG` | Original Phase 3 identifier | Reserved |

#### **Notes:**
- **(31)** Mandatory when ancillary files are provided in a non-FITS format.
- **(33)** Flux maps must always be associated with the RMS noise map, or the SNR map, or both.
- **(37)** If `RADESYS='FK5'`, then `EQUINOX=2000.0` is mandatory.
- **(38)** Mandatory if the system used is other than UTC.
- **(40)** If a refereed publication is not available at the time of the data release, this value can be left empty.
- **(50)** Mandatory if sky resolution is expected to vary within the data collection due to estimation methods.

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
| `ARCFILE` | Archive file name | Reserved |
| `CHECKSUM` | FITS checksum | Mandatory |
| `DATASUM` | FITS data checksum | Mandatory |
| `ORIGFILE` | Original file name | Reserved |
| `P3ORIG` | Original Phase 3 identifier | Reserved |

## 4. File Structure and Validation
- **Simple FITS format**: Science data in primary HDU.
- **Multi-Extension FITS (MEF)**:
  - Science data in separate HDU.
  - Error, confidence, and quality maps stored in additional HDUs.
- Keywords for MEF files:
  - `SCIDATA`: Science data extension name
  - `ERRDATA`: Error data extension name
  - `QUALDATA`: Quality map extension name

## 5. Summary of Best Practices
- Follow the **FITS keyword requirements** for WCS, flux calibration, and data quality.
- Use **Multi-Extension FITS (MEF)** for flux maps with multiple layers.
- Ensure metadata consistency and validation before submission.
- Utilize **ESO-provided validation tools** to check compliance.

For further details and example FITS headers, refer to the [ESO Phase 3 FAQ](https://www.eso.org/sci/observing/phase3/faq.html) and the [ESO Phase 3 Overview](https://www.eso.org/sci/observing/phase3/overview.html).