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

### 3.1 Primary HDU Keywords

These keywords are required in the primary header to capture overall observation and instrument metadata.

| **Keyword**   | **Description**                                                                       | **Location**   | **Requirement**              |
|---------------|---------------------------------------------------------------------------------------|----------------|------------------------------|
| `ORIGIN`      | Data origin institution                                                               | Primary Header | Mandatory                    |
| `TELESCOP`    | Telescope name                                                                        | Primary Header | Mandatory                    |
| `INSTRUME`    | Instrument name                                                                       | Primary Header | Mandatory                    |
| `OBJECT`      | Name of the observed object                                                           | Primary Header | Mandatory                    |
| `DATE-OBS`    | Observation date                                                                      | Primary Header | Mandatory                    |
| `MJD-OBS`     | Modified Julian Date of observation start                                             | Primary Header | Mandatory                    |
| `MJD-END`     | Modified Julian Date of observation end                                               | Primary Header | Mandatory                    |
| `TEXPTIME`    | Total exposure time                                                                   | Primary Header | Mandatory                    |
| `EXPTIME`     | Exposure time per spectral bin                                                        | Primary Header | Mandatory                    |
| `RA`          | Right Ascension of the target                                                         | Primary Header | Mandatory                    |
| `DEC`         | Declination of the target                                                             | Primary Header | Mandatory                    |
| `TIMESYS`     | Time system used for timestamps                                                       | Primary Header | Recommended                  |
| `REFERENC`    | Bibliographic reference identifier                                                    | Primary Header | Recommended                  |
| `CHECKSUM`    | FITS checksum                                                                         | Primary Header | Mandatory                    |
| `DATASUM`     | FITS data checksum                                                                    | Primary Header | Mandatory                    |
| `SPECSYS`     | Spectral reference system                                                             | Primary Header | Mandatory                    |
| `WAVELMIN`    | Minimum wavelength coverage                                                           | Primary Header | Mandatory                    |
| `WAVELMAX`    | Maximum wavelength coverage                                                           | Primary Header | Mandatory                    |
| `PRODCATG`    | Data product category                                                                 | Primary Header | Mandatory                    |
| `RADESYS`     | Reference frame for coordinates (e.g., `FK5`, `ICRS`)                                 | Primary Header | Mandatory                    |
| `EQUINOX`     | Equinox of the coordinate system (e.g., must be 2000.0 if `RADESYS` is `FK5`)           | Primary Header | Mandatory When Applicable    |
| `ASSOCi`      | Association identifier for related data products (e.g., linking ancillary data)         | Primary Header | Mandatory When Applicable    |
| `OBSTECH`     | Observation technique used to acquire the data                                          | Primary Header | Mandatory When Applicable    |
| `EXT_OBJ`     | External object identifier for externally submitted data products                       | Primary Header | Mandatory When Applicable    |

### 3.2 Extension HDU Keywords

These keywords are used in the extension header(s) and provide detailed spectral, calibration, and error information. They are essential for Multi-Extension FITS (MEF) files where additional data layers are provided.

| **Keyword**   | **Description**                                                             | **Location**     | **Requirement** |
|---------------|-----------------------------------------------------------------------------|------------------|-----------------|
| `NELEM`       | Number of elements in the data array                                        | Extension Header | Mandatory       |
| `VOCLASS`     | Virtual Observatory classification of the dataset                         | Extension Header | Mandatory       |
| `VOPUB`       | Virtual Observatory publication status or identifier                       | Extension Header | Mandatory       |
| `TITLE`       | Title or name of the extension data product                                 | Extension Header | Mandatory       |
| `APERTURE`    | Aperture size or setting used during the observation                        | Extension Header | Mandatory       |
| `TELAPSE`     | Elapsed time for the observation (duration)                                 | Extension Header | Mandatory       |
| `TMID`        | Midpoint time of the observation                                            | Extension Header | Mandatory       |
| `SPEC_VAL`    | Central spectral value                                                      | Extension Header | Mandatory       |
| `SPEC_BW`     | Spectral bandwidth                                                          | Extension Header | Mandatory       |

## Additional Considerations

- **Consistency Across HDUs:**  
  Keywords such as `RA`, `DEC`, `EXPTIME`, `CHECKSUM`, and `DATASUM` must have identical values in both the primary and extension headers to maintain consistency.

- **MEF File Structure:**  
  For Multi-Extension FITS (MEF) files, include keywords like `SCIDATA`, `ERRDATA`, and `QUALDATA` when your dataset includes additional layers (e.g., separate error arrays or quality maps).

- **Conditional Keywords:**  
  Keywords marked as “Mandatory When Applicable” must be included based on the specific characteristics of your dataset. Review the Phase 3 documentation for guidance on when these conditions apply.

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

