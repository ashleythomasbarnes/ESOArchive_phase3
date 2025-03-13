**Phase 3 Science Data Products (SDP) Submission Guide for Imaging Products**

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **imaging products**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

Imaging products are science-ready images derived from astronomical observations, including stacked, combined, or processed images. These must follow specific formatting and metadata requirements to be valid for ESO archive ingestion.

## 2. Imaging Product Submission Requirements

### 2.1 General Data Format
- Imaging products must be provided in **FITS format** with appropriate metadata.
- Supported formats:
  - **Single images**: Simple FITS (SCIENCE.IMAGE)
  - **Mosaic images**: Multi-Extension FITS (SCIENCE.MEFIMAGE)
- The **BUNIT** keyword must specify the physical units of the image (e.g., ADU, e-/s).
- All imaging data must include associated error and quality maps.

### 2.2 Astrometry Requirements
- World Coordinate System (WCS) must follow the **FITS convention**.
- **Celestial reference system**: International Celestial Reference System (ICRS).
- **Astrometric error quantification**:
  - `CRDERi`: Random error
  - `CSYERi`: Systematic error
- `RADESYS` must be set to either 'ICRS' or 'FK5'. If 'FK5' is used, `EQUINOX` must be `2000.0`.
- RA/DEC values must be in the correct range: `[0,360]`, `[-90,90]`.

### 2.3 Photometry Requirements
- **Flux scale** should be specified in either:
  - **Johnson system** (normalized to Vega)
  - **AB photometric system**
- **Photometric zero point**: `PHOTZP`
- **Uncertainty of photometric calibration**: `PHOTZPER`
- **Dynamic range characterization**:
  - `ABMAGLIM`: Limiting magnitude
  - `ABMAGSAT`: Saturation magnitude
- `FLUXCAL` must be either 'ABSOLUTE' or 'UNCALIBRATED'.

### 2.4 Temporal Information
- Observation time must be recorded using **Modified Julian Date (MJD)**.
- `MJD-END` must be greater than or equal to `MJD-OBS`.
- `EXPTIME` (exposure time per pixel) and `TEXPTIME` (total exposure time) must be included.

## 3. Required and Recommended FITS Header Keywords
The following keywords must be present in the FITS headers for imaging data:

### 3.1 General Information
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `ORIGIN` | Data origin institution. | Primary Header | Mandatory |
| `TELESCOP` | Telescope name. | Primary Header | Mandatory |
| `INSTRUME` | Instrument name. | Primary Header | Mandatory |
| `OBJECT` | Name of the observed object. | Primary Header | Mandatory |
| `FILTER` | Filter used during observation. | Primary Header | Mandatory |
| `DATE-OBS` | Observation date. | Primary Header | Mandatory |
| `MJD-OBS` | Modified Julian Date of observation start. | Primary Header | Mandatory |
| `MJD-END` | Modified Julian Date of observation end. | Primary Header | Mandatory |
| `EXPTIME` | Total exposure time per pixel. | Primary & Extension Header | Mandatory |
| `TEXPTIME` | Total exposure time. | Primary Header | Mandatory |
| `RA` | Right Ascension of the target. | Primary & Extension Header | Mandatory |
| `DEC` | Declination of the target. | Primary & Extension Header | Mandatory |
| `TIMESYS` | Time system used for timestamps. | Primary Header | Recommended |
| `REFERENC` | Bibliographic reference identifier. | Primary Header | Recommended |
| `NAXIS1`, `NAXIS2` | Image dimensions. | Primary Header | Mandatory |
| `CHECKSUM` | FITS checksum. | Primary & Extension Header | Mandatory |
| `DATASUM` | FITS data checksum. | Primary & Extension Header | Mandatory |

### 3.2 WCS and Astrometry
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `RADESYS` | Coordinate reference system (e.g., 'ICRS'). | Primary & Extension Header | Mandatory |
| `CRVALi` | Coordinate value at reference pixel. | Extension Header | Mandatory |
| `CRPIXi` | Reference pixel coordinates. | Extension Header | Mandatory |
| `CTYPEi` | Coordinate system type (e.g., 'RA---TAN', 'DEC--TAN'). | Extension Header | Mandatory |
| `CUNITi` | Units of coordinates (typically 'deg'). | Extension Header | Mandatory |
| `CDi_j`  | WCS transformation matrix elements. | Extension Header | Mandatory |

### 3.3 Photometry
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `FLUXCAL` | Flux calibration status ('ABSOLUTE' or 'UNCALIBRATED'). | Extension Header | Mandatory |
| `BUNIT` | Physical unit of flux (e.g., ADU, e-/s). | Extension Header | Mandatory |
| `PHOTZP`  | Photometric zero point. | Extension Header | Recommended |
| `PHOTZPER` | Photometric calibration uncertainty. | Extension Header | Recommended |
| `ABMAGLIM` | Limiting magnitude. | Extension Header | Recommended |
| `ABMAGSAT` | Saturation magnitude. | Extension Header | Recommended |

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
- Follow the **FITS keyword requirements** for WCS, photometry, and data quality.
- Use **Multi-Extension FITS (MEF)** for imaging data with multiple layers.
- Ensure metadata consistency and validation before submission.
- Utilize **ESO-provided validation tools** to check compliance.

For further details and example FITS headers, refer to the [ESO Phase 3 FAQ](https://www.eso.org/sci/observing/phase3/faq.html) and the [ESO Phase 3 Overview](https://www.eso.org/sci/observing/phase3/overview.html).

