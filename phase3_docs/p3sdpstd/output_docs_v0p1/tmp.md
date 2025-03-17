**Phase 3 Science Data Products (SDP) Submission Guide for Spectral Data**

## 1. Introduction
The Phase 3 process facilitates the submission, validation, and ingestion of science data products (SDPs) into the ESO Science Archive. This document provides a streamlined guide for submitting **spectral data**, ensuring compliance with ESO/SDP standards while simplifying the process for users.

Spectral data products include one-dimensional spectra, extracted spectra from integral field units, and multi-object spectroscopy. These must follow strict formatting and metadata requirements to ensure archive compatibility.

## 3. Generic FITS Standard Keywords

### 3.1 Primary HDU Keywords (Generic FITS Standard)
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

### 3.2 Extension HDU Keywords (Generic FITS Standard)
| Keyword  | Description | Type |
|----------|-------------|------|
| `XTENSION` | Type of FITS extension (`IMAGE`, `TABLE`, `BINTABLE`) | Mandatory |
| `BITPIX`   | Number of bits per pixel in the extension | Mandatory |
| `NAXIS`    | Number of data axes in the extension | Mandatory |
| `NAXIS1`   | Length of first data axis in the extension | Mandatory |
| `NAXIS2`   | Length of second data axis in the extension | Mandatory |
| `PCOUNT`   | Parameter count (used in binary tables) | Mandatory |
| `GCOUNT`   | Group count (typically `1` for standard FITS files) | Mandatory |
| `EXTNAME`  | Name of the extension HDU | Recommended |
| `CHECKSUM` | FITS checksum to verify data integrity | Mandatory |
| `DATASUM`  | Data checksum to verify extension contents | Mandatory |

These keywords form the **essential structure of any FITS file** and are required by the FITS standard. Let me know if you need any refinements! 🚀

