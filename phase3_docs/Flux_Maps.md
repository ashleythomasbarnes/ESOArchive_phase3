**Phase 3 Science Data Products (SDP) Submission Guide for Flux Maps**

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **flux maps**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

Flux maps are specialized imaging products where pixel values represent flux measurements in physical units, often derived from calibrated observations. They must follow specific formatting and metadata requirements to be valid for ESO archive ingestion.

## 2. Flux Map Submission Requirements

### 2.1 General Data Format
- Flux maps must be provided in **FITS format** with appropriate metadata.
- Supported formats:
  - **Single images**: Simple FITS (SCIENCE.IMAGE.FLUXMAP)
  - **Multi-extension images**: Multi-Extension FITS (MEF)
- The **BUNIT** keyword must specify the physical units of flux (e.g., Jy/pixel, W/m²/sr, or MJy/sr).
- All flux maps must include associated error and quality maps.

### 2.2 Astrometry Requirements
- World Coordinate System (WCS) must follow the **FITS convention**.
- **Celestial reference system**: International Celestial Reference System (ICRS).
- **Astrometric error quantification**:
  - `CRDERi`: Random error
  - `CSYERi`: Systematic error
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
The following keywords must be present in the FITS headers for flux map data:

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
| `MAPMODE` | Mapping mode used for observation. | Primary Header | Recommended |
| `SKY_RES` | Sky resolution. | Primary Header | Recommended |
| `BNOISE` | Background noise level. | Primary Header | Recommended |
| `NAXIS1`, `NAXIS2` | Image dimensions. | Primary Header | Mandatory |
| `WAVELMIN`, `WAVELMAX` | Minimum and maximum wavelengths. | Primary Header | Mandatory |
| `NCOMBINE` | Number of combined exposures. | Primary Header | Recommended |
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

### 3.3 Flux Calibration and Photometry
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `FLUXCAL` | Flux calibration status ('ABSOLUTE' or 'UNCALIBRATED'). | Extension Header | Mandatory |
| `BUNIT` | Physical unit of flux (e.g., Jy/pixel, MJy/sr). | Extension Header | Mandatory |
| `FLUXERR` | Flux error per pixel. | Extension Header | Recommended |
| `PHOTZP`  | Photometric zero point. | Extension Header | Recommended |

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

