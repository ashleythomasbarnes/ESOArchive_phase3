**Phase 3 Science Data Products (SDP) Submission Guide for Integral Field Spectroscopy (IFS) Data Cubes**

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **IFS data cubes**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

IFS data cubes contain three-dimensional spectral data where two axes correspond to spatial coordinates, and the third represents wavelength. These products require specific metadata and formatting to be valid for ESO archive ingestion.

## 2. IFS Cube Submission Requirements

### 2.1 General Data Format
- IFS cubes must be provided in **FITS format** with appropriate metadata.
- Supported formats:
  - **Integral Field Spectroscopy Cube**: SCIENCE.CUBE.IFS
- The **BUNIT** keyword must specify the physical units of flux (e.g., erg/s/cm²/Å, W/m²/Hz).
- All data cubes must include associated error, variance, and quality maps.

### 2.2 Spatial and Wavelength Calibration
- World Coordinate System (WCS) must follow the **FITS convention**.
- **Celestial reference system**: International Celestial Reference System (ICRS).
- **Astrometric accuracy**:
  - `CRDERi`: Random error
  - `CSYERi`: Systematic error
- **Wavelength Calibration**:
  - Wavelength axis must be defined in **Angstrom (Å), nanometers (nm), or microns (µm)**.
  - `SPEC_BIN`: Spectral bin width.
  - `SPEC_RES`: Spectral resolution.
  - `WAVELMIN`, `WAVELMAX`: Minimum and maximum wavelengths covered.

### 2.3 Flux Calibration and Photometry Requirements
- **Flux calibration status** must be specified using `FLUXCAL`:
  - 'ABSOLUTE' if the cube is fully flux-calibrated.
  - 'UNCALIBRATED' if only relative flux calibration is applied.
- **Flux uncertainty quantification**:
  - `FLUXERR`: Error associated with each spectral bin.
  - `SPEC_ERR`: Spectral uncertainty.

### 2.4 Temporal Information
- Observation time must be recorded using **Modified Julian Date (MJD)**.
- `MJD-END` must be greater than or equal to `MJD-OBS`.
- `EXPTIME` (exposure time per spectral bin) and `TEXPTIME` (total exposure time) must be included.

## 3. Required and Recommended FITS Header Keywords
The following keywords must be present in the FITS headers for IFS data cubes:

### 3.1 General Information
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `ORIGIN` | Data origin institution. | Primary Header | Mandatory |
| `TELESCOP` | Telescope name. | Primary Header | Mandatory |
| `INSTRUME` | Instrument name. | Primary Header | Mandatory |
| `OBJECT` | Name of the observed object. | Primary Header | Mandatory |
| `DATE-OBS` | Observation date. | Primary Header | Mandatory |
| `MJD-OBS` | Modified Julian Date of observation start. | Primary Header | Mandatory |
| `MJD-END` | Modified Julian Date of observation end. | Primary Header | Mandatory |
| `EXPTIME` | Total exposure time per bin. | Primary & Extension Header | Mandatory |
| `TEXPTIME` | Total exposure time. | Primary Header | Mandatory |
| `RA` | Right Ascension of the target. | Primary & Extension Header | Mandatory |
| `DEC` | Declination of the target. | Primary & Extension Header | Mandatory |
| `TIMESYS` | Time system used for timestamps. | Primary Header | Recommended |
| `REFERENC` | Bibliographic reference identifier. | Primary Header | Recommended |
| `CHECKSUM` | FITS checksum. | Primary & Extension Header | Mandatory |
| `DATASUM` | FITS data checksum. | Primary & Extension Header | Mandatory |

### 3.2 WCS and Astrometry
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `RADESYS` | Coordinate reference system (e.g., 'ICRS'). | Primary & Extension Header | Mandatory |
| `CRVALi` | Coordinate value at reference pixel. | Extension Header | Mandatory |
| `CRPIXi` | Reference pixel coordinates. | Extension Header | Mandatory |
| `CTYPEi` | Coordinate system type (e.g., 'RA---TAN', 'DEC--TAN', 'WAVE'). | Extension Header | Mandatory |
| `CUNITi` | Units of coordinates (typically 'deg' for spatial and 'Angstrom' for spectral). | Extension Header | Mandatory |
| `CDi_j`  | WCS transformation matrix elements. | Extension Header | Mandatory |

### 3.3 Spectral and Flux Calibration
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `SPECSYS` | Spectral reference system. | Primary Header | Mandatory |
| `WAVELMIN`, `WAVELMAX` | Minimum and maximum wavelength coverage. | Primary Header | Mandatory |
| `SPEC_RES` | Spectral resolution. | Extension Header | Mandatory |
| `SPEC_BIN` | Spectral bin width. | Extension Header | Recommended |
| `FLUXCAL` | Flux calibration status ('ABSOLUTE' or 'UNCALIBRATED'). | Extension Header | Mandatory |
| `BUNIT` | Physical unit of flux (e.g., erg/s/cm²/Å). | Extension Header | Mandatory |
| `FLUXERR` | Flux error per bin. | Extension Header | Recommended |
| `SPEC_ERR` | Spectral uncertainty. | Extension Header | Recommended |

## 4. File Structure and Validation
- **FITS Format for Data Cubes**:
  - Three-dimensional array (spatial-spatial-spectral).
  - First two axes represent celestial coordinates, third axis represents wavelength.
- **Multi-Extension FITS (MEF) Format**:
  - Science data in separate HDU.
  - Error estimates, variance, and quality maps stored in additional HDUs.
- Keywords for MEF files:
  - `SCIDATA`: Science data extension name
  - `ERRDATA`: Error data extension name
  - `QUALDATA`: Quality map extension name

## 5. Summary of Best Practices
- Follow the **FITS keyword requirements** for WCS, spectral calibration, flux calibration, and data quality.
- Use **Multi-Extension FITS (MEF)** for IFS data cubes with multiple layers.
- Ensure metadata consistency and validation before submission.
- Utilize **ESO-provided validation tools** to check compliance.

For further details and example FITS headers, refer to the [ESO Phase 3 FAQ](https://www.eso.org/sci/observing/phase3/faq.html) and the [ESO Phase 3 Overview](https://www.eso.org/sci/observing/phase3/overview.html).

