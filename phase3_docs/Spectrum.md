**Phase 3 Science Data Products (SDP) Submission Guide for Spectral Data**

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **spectral data**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

Spectral data products include one-dimensional spectra, extracted spectra from integral field units, and multi-object spectroscopy. These must follow strict formatting and metadata requirements to ensure archive compatibility.

## 2. Spectral Data Submission Requirements

### 2.1 General Data Format
- Spectral data must be provided in **FITS format** with appropriate metadata.
- Supported formats:
  - **Single spectrum**: SCIENCE.SPECTRUM (1D extracted spectrum)
  - **Multi-object spectrum**: SCIENCE.MOSSPECTRA
  - **Integral Field Unit (IFU) spectrum**: SCIENCE.CUBE.IFS
- The **BUNIT** keyword must specify the physical units of flux (e.g., erg/s/cm²/Å, W/m²/Hz).
- All spectral data must include associated error arrays, noise estimates, and quality flags where applicable.

### 2.2 Wavelength Calibration
- Spectral axes must be correctly defined using **FITS WCS conventions for spectral data**.
- **Coordinate system for wavelengths**:
  - Wavelength must be given in **Angstrom (Å), nanometers (nm), or microns (µm)**.
  - Frequency units (GHz) or energy (keV) are acceptable for certain datasets.
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

### 2.4 Temporal Information
- Observation time must be recorded using **Modified Julian Date (MJD)**.
- `MJD-END` must be greater than or equal to `MJD-OBS`.
- `EXPTIME` (exposure time per spectral bin) and `TEXPTIME` (total exposure time) must be included.

## 3. Required and Recommended FITS Header Keywords
The following keywords must be present in the FITS headers for spectral data:

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

### 3.2 Wavelength Calibration and Spectral Characteristics
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `SPECSYS` | Spectral reference system. | Primary Header | Mandatory |
| `SPEC_VAL` | Central spectral value. | Extension Header | Mandatory |
| `SPEC_BW` | Spectral bandwidth. | Extension Header | Mandatory |
| `WAVELMIN`, `WAVELMAX` | Minimum and maximum wavelength coverage. | Primary Header | Mandatory |
| `SPEC_RES` | Spectral resolution. | Extension Header | Mandatory |
| `SPEC_BIN` | Spectral bin width. | Extension Header | Recommended |

### 3.3 Flux Calibration and Uncertainty
| Keyword  | Description | Location | Type |
|----------|-------------|----------|------|
| `FLUXCAL` | Flux calibration status ('ABSOLUTE' or 'UNCALIBRATED'). | Extension Header | Mandatory |
| `BUNIT` | Physical unit of flux (e.g., erg/s/cm²/Å). | Extension Header | Mandatory |
| `FLUXERR` | Flux error per bin. | Extension Header | Recommended |
| `SPEC_ERR` | Spectral uncertainty. | Extension Header | Recommended |

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

