**Spectrum** - Phase 3 Science Data Products (SDP) Submission Guide

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **spectral data**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

Spectral data products include one-dimensional spectra from single-object and multi-object spectroscopic observations. These must follow strict formatting and metadata requirements to ensure archive compatibility.

## 3. Required and Recommended FITS Header Keywords

### 3.2.2 Extension HDU Keywords

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

### 3.2.3 Extension HDU Table Column Definitions

| `TTYPEi` Value | Description |
|---------------|-------------|
| `WAVE` | The wavelength array |
| `FREQ` | The frequency array |
| `ENER` | The energy array |
| `FLUX` | The data spectrum: either the sky-background subtracted spectrum or the continuum normalized spectrum. |
| `ERR` | The error spectrum. Errors must be provided in the same units as the flux array and cannot be expressed as a percentage. |
| `QUAL` | An array of integer values: 0 = good data, 1 = bad (unspecified reason), other positive integers flag bad or dubious data. If absent, all values are assumed good. Encoding should use powers of 2 for quality conditions. |
| `BGFLUX` | The sky background spectrum |
| `CONTINUUM` | The continuum spectrum |
| `EXPOSURE` | The exposure array for combined spectra of different wavelength coverage |

These definitions ensure proper structuring and compliance with ESO Phase 3 standards for spectral data submissions.

