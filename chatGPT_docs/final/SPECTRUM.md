# **Spectrum** - Phase 3 Science Data Products (SDP) Submission Guide

## 1. Introduction

**Purpose:** This guide condenses the ESO SDP standard and DICD into clear, practical steps for preparing and submitting 1D spectra to Phase 3. The **data category** you must write in the header for a spectrum is:

```
PRODCATG = 'SCIENCE.SPECTRUM'
```

**What counts as SPECTRUM?** A **single-target 1D spectrum** delivered as a **single-row FITS BINTABLE** with all required metadata and consistent arrays.

**Key principles.**

* Products are fully calibrated, documented, and traceable.
* FITS headers capture observation context, calibration status, and VO interoperability.

### 1.1 What is Phase 3?

**Phase 3** is ESO’s process to **prepare, validate, ingest, and publish** Science Data Products (SDPs) in the ESO Science Archive. SDPs are **fully calibrated** (instrumental/atmospheric signatures removed), in **physical units**, with **documented noise properties** (e.g. S/N, limiting magnitude).

**Who produces SDPs:**
- PIs of ESO programmes (public surveys, large and calibration programmes)
- ESO scientists via pipelines/QC or dedicated re-processing of homogeneous data sets
- Users of **non-ESO telescopes** (e.g. GTC, NGTS)
- PIs of normal ESO programmes (voluntary)
- Community members publishing **archival** results

**Policy:** Phase 3 is **mandatory** for **Public Surveys** and **Large Programmes** since **Period 75**; **optional** for other programmes.

**Support:** ESO provides **standards**, **procedures/infrastructure** (Release Manager, FTP), and **tools** for data preparation.

**Audience:** (1) PIs/collaborators returning reduced data from ESO or non-ESO facilities; (2) ESO scientists in QC/reprocessing; (3) instrument scientists & pipeline developers; (4) archive users who need to understand SDP structure/format.

### 1.2 Overview of how to submit data to Phase 3

**1) Register your collection.** Open the **Phase 3 Release Manager (RM)** and log in with **ESO User Portal** credentials. If your programme isn’t listed under *Data Collections*, contact Phase 3 Operations via the Helpdesk with subject `REQUEST FOR PHASE 3 PROGRAMME <PPP.C-NNNN>`. You’ll receive confirmation (typically within one working day) and the collection will appear in RM.

**2) Prepare your data.** Format files per the **ESO Science Data Products Standard (latest)** -- see Section 2 below for SPECTRUM specifics.

**3) Upload your data.** FTP to `phase3ftp.eso.org` using the **path shown in RM**: `/<Data collection name>/batch_<ID>`. Use any client (e.g. **lftp**, **FileZilla**); log in with your User Portal credentials. Access is limited to the **PI**, **survey manager**, and **delegates**. Files can be replaced/removed **until you close** the batch. Resolve format issues first (see automatic **Phase 3 checks**).

**4) Run server-side checks.** In RM, click **CLOSE** to trigger format/provenance verification. Processing time depends on volume; you’ll get an **e‑mail** when done. Fix any issues and re-open if needed. *Note:* closing makes the FTP directory **read‑only**.

**5) Upload the release description (PDF).** Summarise **release content**, **originating observations**, **calibration & reduction**, **data quality**, **data format**, and (optionally) **scientific context**. For more information on preparing this document, see the [Release Description Overview](https://www.eso.org/sci/observing/phase3/release_description.html) page, and for the template, see the [Phase 3 Release Description Template](https://www.eso.org/sci/observing/phase3/release-description-tmpl.docx) page. 

**6) Finalise and submit.** The **PI** (not delegable) presses **SUBMIT** to confirm consistency with the release description. Review the **Content Summary** (file counts, data types, dates, sky coverage). ESO performs **final content validation** and publishes the data; if issues remain, you’ll receive a **report** with next steps.

---

## 2. Spectral Data Submission Requirements

### 2.1 General Data Format

- Supported format:
  - **SCIENCE.SPECTRUM**: 1D extracted spectrum.
- Spectral data must be provided in **FITS format** with appropriate metadata. 
  - A **Primary Header** with no data (`NAXIS=0`).
  - A **single Extension** containing a **BINTABLE** with `NAXIS=2`, where the data arrays are stored as vectors in single cells, meaning `NAXIS2=1`.
  - Each field of the **BINTABLE** must be further described in the extension header.
  - Only **standard BINTABLE extensions** are allowed, implying `PCOUNT=0` and `GCOUNT=1`.
- Information associated with the science spectrum shall be stored within the **same extension** as the main science spectrum, including:
  - Spectrum - `FLUX`
  - Error spectrum - `ERR`
  - Data quality information - `QUAL`
  - Sky background - `BGFLUX`
  - Best-fitted model for the continuum - `CONTINUUM`
  - The exposure array for e.g. combined spectra of different wavelength
coverage - `EXPOSURE`
- 2D spectral frames may be submitted in addition as associated files.

### 2.2 Spectral Data Model and VO Compliance

- **Case 1: Single spectral coordinate array and single flux array**
  - The **IVOA Spectrum Data Model v1.1** must be used to describe arrays.
  - The `VOCLASS` keyword in the first extension header must be set to `'SPECTRUM v1.0'`.
  - The following keywords must be defined:
    ```
    VOCLASS = 'SPECTRUM v1.0' / VO Data Model
    TTYPE1  = 'WAVE '  / Label for field 1
    TTYPE2  = 'FLUX '  / Label for field 2
    TTYPE3  = 'ERR '   / Label for field 3
    TUTYP1  = 'Spectrum.Data.SpectralAxis.Value'
    TUTYP2  = 'Spectrum.Data.FluxAxis.Value'
    TUTYP3  = 'Spectrum.Data.FluxAxis.Accuracy.StatError'
    ```

- **Case 2: Multiple flux arrays computed with different recipes**
  - The **IVOA Spectrum Data Model v2.0** must be used.
  - The `VOCLASS` keyword must be set to `'SPECTRUM v2.0'`.
  - The submitter must select which flux column is the **default flux**.
  - Example:
    ```
    VOCLASS = 'SPECTRUM v2.0' / VO Data Model
    TTYPE1  = 'WAVE '  / Label for field 1
    TTYPE2  = 'FLUX_A '  / Flux from recipe A
    TTYPE3  = 'ERR_A '   / Error for FLUX_A
    TTYPE4  = 'FLUX_B '  / Flux from recipe B
    TTYPE5  = 'ERR_B '   / Error for FLUX_B
    TUTYP1  = 'spec:Data.SpectralAxis.Value'
    TUTYP2  = 'spec:Data.FluxAxis.Value'
    TUTYP3  = 'spec:Data.FluxAxis.Accuracy.StatError'
    TUTYP4  = 'eso:Data.FluxAxis.Value'
    TUTYP5  = 'eso:Data.FluxAxis.Accuracy.StatError'
    ```

- **Continuum Normalized Spectra**
  - If a normalized spectrum is published as the main flux, but the unnormalized flux is also stored, the **CONTNORM** keyword must be included in the **Primary Header**:
    ```
    CONTNORM = T
    ```
  - Example extension header:
    ```
    VOCLASS = 'SPECTRUM v2.0' / VO Data Model
    TTYPE1 = 'WAVE ' / Label for field 1
    TTYPE2 = 'FLUX_NORM' / Label for field 2
    TTYPE3 = 'ERR_NORM' / Label for field 3
    TTYPE4 = 'FLUX ' / Label for field 4
    TTYPE5 = 'ERR ' / Label for field 5
    TUTYP1 = 'spec:Data.SpectralAxis.Value'
    TUTYP2 = 'spec:Data.FluxAxis.Value'
    TUTYP3 = 'spec:Data.FluxAxis.Accuracy.StatError'
    TUTYP4 = 'eso:Data.FluxAxis.Value'
    TUTYP5 = 'eso:Data.FluxAxis.Accuracy.StatError'
    TUCD1 = 'em.wl;obs.atmos' / em.wl: vacuum, em.wl;obs.atmos: air
    TUCD2 = 'phot.flux.density;em.wl;arith.ratio;meta.main'
    TUCD3 = 'stat.error;phot.flux.density;em.ql;arith.ratio;meta.main'
    TUCD4 = 'phot.flux.density;em.wl'
    TUCD5 = 'stat.error;phot.flux.density;em.ql'
    ```

### 2.3 Additional Requirements
- All data arrays in the **first row of the binary table must have the same number of points** (`NELEM`).
- All `EXTNAME` keyword values **must be unique** within a given FITS file.
- If present, there may be at most one `BGFLUX` field.
- For data not normalized to the continuum, `TUCD` and `TUNIT` of `BGFLUX` must match those of `FLUX`.
- To differentiate between air and vacuum wavelengths, `TUCD1` must be set accordingly:
  - **Air Wavelength:**
    ```
    TTYPE1 = 'WAVE'
    TUCD1  = 'em.wl;obs.atmos' / Air wavelength
    ```
  - **Vacuum Wavelength:**
    ```
    TTYPE1 = 'WAVE'
    TUCD1  = 'em.wl' / Vacuum wavelength
    ```
  - `WAVELMIN` and `WAVELMAX` must be consistent with the spectral axis definition.

---

## 3. Required and Recommended FITS Header Keywords

#### 3.1.1 Primary HDU Keywords (Generic FITS Standard)

| Keyword  | Description | Type |
|----------|-------------|------|
| `SIMPLE`  | Conforms to FITS standard (must be `T`) | Mandatory |
| `BITPIX`  | Number of bits per pixel (e.g., `-32` for floating-point data) | Mandatory |
| `NAXIS`   | Number of data axes | Mandatory |
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
| `ORIGIN` | Observatory or facility where the data were originally obtained | Mandatory |  
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
| `EXT_OBJ` | Extended object (TRUE), Pointlike (FALSE) | Mandatory (42) |  
| `CONTNORM` | Continuum normalization factor | Mandatory |  
| `TOT_FLUX` | Total flux flag (TRUE/FALSE) | Mandatory |  
| `FLUXERR` | Global percentage scale error | Mandatory (43) |  
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
| `NOESODAT` | Non-ESO proprietary data flag | Mandatory (59) |  

#### **Notes:**

- **(31)** Mandatory when ancillary files are provided in a non-FITS format.
- **(32)** Mandatory when ancillary files are associated with the scientific data.
- **(37)** If RADESYS='FK5', then EQUINOX=2000.0 is mandatory.
- **(38)** Mandatory if the system used is other than UTC.
- **(40)** If a refereed publication is not available at the time of the data release, this value can be left empty.
- **(42)** Mandatory for externally submitted data products.
- **(43)** Applies to `FLUXCAL='ABSOLUTE'` (must be case if `TOT_FLUX=TRUE`). If `FLUXCAL='UNCALIBRATED'`, set `FLUXERR=-1`.
- **(51)** Mandatory for adaptive optics (AO) observations.
- **(59)** Set NOESODAT=T for non-ESO proprietary data.

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
| `OBJECT` | Target name (must match primary) | Mandatory |   
| `RA`, `DEC` | Right Ascension and Declination | Mandatory |  
| `NELEM` | Number of elements in the data array | Mandatory |  
| `VOCLASS` | Virtual Observatory classification | Mandatory |  
| `VOPUB` | Name of the publisher - `'ESO/SAF'` | Mandatory |  
| `TITLE` | Title of the dataset | Mandatory |  
| `APERTURE` | Aperture size used during observation in degrees | Mandatory |  
| `TELAPSE` | Elapsed observation time | Mandatory |  
| `TMID` | Midpoint time of observation | Mandatory |  
| `SPEC_VAL` | Central spectral value | Mandatory |  
| `SPEC_BW` | Spectral bandwidth | Mandatory |  
| `ARCFILE` | Archive file name | Reserved |  
| `CHECKSUM` | FITS checksum | Mandatory |  
| `DATASUM` | FITS data checksum | Mandatory |  
| `ORIGFILE` | Original file name | Reserved |  
| `P3ORIG` | Original Phase 3 identifier | Reserved |  
| `TFIELDS` | Number of columns in binary table | Mandatory |  
| `TTYPEi` | Column name | Mandatory |  
| `TFORMi` | Data format of column | Mandatory |  
| `TCOMMi` | Comment on column usage | Optional |  
| `TUNITi` | Units of column values | Mandatory (61) |  
| `TUTYPi` | Utype of column | Mandatory (62) |  
| `TUCDi` | Unified Content Descriptor (UCD) | Mandatory |  
| `TDMINi` | Minimum value in data array | Mandatory (63) |  
| `TDMAXi` | Maximum value in data array | Mandatory (63) |  
| `EXTNAME` | Extension name | Mandatory (65) |  

#### **Notes:**

- **(61)** The unit of each column must follow SI standards and be properly defined using TUNITi.
- **(62)** TUTYPi must be set following IVOA Utype conventions (e.g., 'Spectrum.Data.SpectralAxis.Value').
- **(63)** TDMINi and TDMAXi should be specified when applicable to define valid data ranges in the column.
- **(65)** EXTNAME must be unique within a given FITS file.

#### 3.2.2 `TTYPEi` keyword values 

Example usage in header: `TTYPE1 = 'WAVE'` for spectral axis given in wavelength.

| `TTYPEi` Value | Description |  
|---------------|-------------|  
| `WAVE` | The wavelength array |  
| `FREQ` | The frequency array |  
| `ENER` | The energy array |  
| `FLUX` | The data spectrum: either the sky-background subtracted spectrum or the continuum-normalized spectrum. |  
| `ERR` | The error spectrum. Errors must be provided in the same units as the flux array and cannot be expressed as a percentage. |  
| `QUAL` | An array of integer values: 0 = good data, 1 = bad (unspecified reason), other positive integers flag bad or dubious data. If absent, all values are assumed good. Encoding should use powers of 2 for quality conditions. |  
| `BGFLUX` | The sky background spectrum |  
| `CONTINUUM` | The continuum spectrum |  
| `EXPOSURE` | The exposure array for combined spectra of different wavelength coverage | 


#### Notes for Extension HDU Table Column Definitions  

- **Wavelength Units:**  
  - When measured in air: `TUCD1 = 'em.wl;obs.atmos'`  
  - When measured in vacuum: `TUCD1 = 'em.wl'`  
  - `WAVELMIN` and `WAVELMAX` must be consistent with the spectral axis definition.  

- **Error Representation:**  
  - `ERR` values **must be in the same units as `FLUX`**.  
  - Percentage-based errors **are not allowed**.  

- **Quality Encoding (`QUAL`):**  
  - `0` = Good data  
  - `1` = Bad for an unspecified reason  
  - Values **greater than 1** can be used for specific flags (preferably using powers of 2).  

- **Flux Representation:**  
  - If a **continuum-normalized spectrum** is used, the **CONTNORM** keyword must be set to `T` in the **Primary HDU**.  
  - The **normalized flux** must be stored in the **FLUX_NORM** column, with an additional **FLUX** column for the unnormalized data.  

### Additional Considerations

- **Consistency Across HDUs:**  
  Keywords such as `RA`, `DEC`, `CHECKSUM`, and `DATASUM` must have identical values in both the primary and extension headers to maintain consistency.

- **MEF File Structure:**  
  For Multi-Extension FITS (MEF) files, include keywords like `SCIDATA`, `ERRDATA`, and `QUALDATA` when your dataset includes additional layers (e.g., separate error arrays or quality maps).

- **Validation:**  
  Prior to submission, run your FITS files through the ESO-provided validation tools to ensure that all header keywords comply with Phase 3 SDP standards.

For further details and examples of FITS headers, please refer to the [ESO Phase 3 FAQ](https://www.eso.org/sci/observing/phase3/faq.html) and the [ESO Phase 3 Overview](https://www.eso.org/sci/observing/phase3/overview.html).

## 4. Minimal header skeleton (illustrative)

```fits
SIMPLE  =                T / file does conform to FITS standard
BITPIX  =              -32 / number of bits per data pixel
NAXIS   =                0 / number of data axes
EXTEND  =                T / FITS dataset may contain extensions
COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
DATE    = '2025-10-23T03:26:04.139' / UT date when this file was written
INSTRUME= 'HARPS'          / Instrument used.
RA      =         3.460039 / 00:13:50.4 RA (J2000) pointing (deg)
DEC     =        -52.57105 / -52:34:15.7 DEC (J2000) pointing (deg)
EQUINOX =            2000. / Standard FK5 (years)
RADECSYS= 'FK5'            / Coordinate reference frame
EXPTIME =         299.9959 / Total integration time
MJD-OBS =   60971.13961105 / MJD start (2025-10-23T03:21:02.395)
DATE-OBS= '2025-10-23T03:21:02.395' / Date of observation
UTC     =        12057.000 / 03:20:57.000 UTC at start (sec)
LST     =         2713.705 / 00:45:13.705 LST at start (sec)
PI-COI  = 'UNKNOWN'        / PI-COI name.
OBSERVER= 'UNKNOWN'        / Name of observer.
ORIGIN  = 'ESO'            / European Southern Observatory
TELESCOP= 'ESO-3.6'        / ESO Telescope Name
M_EPOCH =                F / T if resulting from multiple epochs
SINGLEXP=                T / T if resulting from single exposure
PRODLVL =                2 / Product level: 1-raw 2-science grade 3-advanced
DISPELEM= 'HARPS Echelle'  / Dispersive element name
SPECSYS = 'BARYCENT'       / Reference frame
OBJECT  = 'J0013-52'       /   Object
TEXPTIME=         299.9959 / [s] Total integration time of all exposures
EXT_OBJ =                F / T if extended
MJD-END =   60971.14308322 / [d] End of observation
PROG_ID = '116.28LH.001'   / ESO Programme ID
OBID1   =       1000593703 / Observation block ID
PROCSOFT= 'HARPS_3.8'      /         
OBSTECH = 'ECHELLE'        / Observation technique
PRODCATG= 'SCIENCE.SPECTRUM' / Data product category
FLUXCAL = 'UNCALIBRATED'   / Type of flux calibration
CONTNORM=                F / T if normalised to the continuum
WAVELMIN=          378.293 / [nm] Minimum wavelength
WAVELMAX=          691.242 / [nm] Maximum wavelength
SPEC_BIN=            0.001 / [nm] Wavelength bin size
TOT_FLUX=                F / T if photom. cond. and all src flux captured
FLUXERR =             -1.0 / Uncertainty in flux scale (%)
NCOMBINE=                1 / # of combined raw science data files
REFERENC= ' '              / Bibliographic reference
SNR     =             46.1 / Median signal to noise
SPEC_RES=         115000.0 / Reference spectral resolving power
GAIN    =             1.36 / Conversion from electrons to ADU
DETRON  = 4.67247428449605 / Readout noise per output (e-)
ASSOC1  = 'ANCILLARY.HARPSTAR' / Associated file 1 category
ASSON1  = 'HARPS.2025-10-23T03:21:02.395_DRS_HARPS_3.8.tar' / Associated file 1
ASSOM1  = 'fd62d8cfb9003cd4f130b9ca294ab112' / Associated file 1 md5sum
PROV1   = 'HARPS.2025-10-23T03:21:02.395.fits' / Originating file name
CHECKSUM= 'e7S4g4P4e4P4e4P4' / HDU checksum updated 2025-10-24T01:01:11
DATASUM = '0'              / data unit checksum updated 2025-10-24T03:00:32
ARCFILE = 'ADP.2025-10-24T01:01:07.534.fits' / Archive File Name
P3ORIG  = 'IDP'            / ESO internal data product
END
XTENSION= 'BINTABLE'       / binary table extension
BITPIX  =                8 / array data type
NAXIS   =                2 / number of array dimensions
NAXIS1  =          5007200 / length of dimension 1
NAXIS2  =                1 / length of dimension 2
PCOUNT  =                0 / number of group parameters
GCOUNT  =                1 / number of groups
TFIELDS =                3 / number of table fields
TTYPE1  = 'WAVE'           /         
TFORM1  = '312950D'        /         
TTYPE2  = 'FLUX'           /         
TFORM2  = '312950E'        /         
TTYPE3  = 'ERR'            /         
TFORM3  = '312950E'        /         
VOCLASS = 'SPECTRUM v1.0'  / VO Data Model
VOPUB   = 'ESO/SAF'        / VO Publishing Autority
EXTNAME = 'SPECTRUM'       / Extension Name
INHERIT =                T / Primary HDU inherited
TITLE   = 'J0013-52_HARPS.2025-10-23T03:21:02.395_s1d_A' /         
OBJECT  = 'J0013-52'       /   Object
RA      =         3.460039 / 00:13:50.4 RA (J2000) pointing (deg)
DEC     =       -52.571050 / -52:34:15.7 DEC (J2000) pointing (deg)
APERTURE=       0.00027778 / [deg] Aperture diameter
TELAPSE =         299.9959 / [s] Total elapsed time
TMID    =   60971.14134714 / [d] MJD mid exposure
NELEM   =           312950 / Length of the data arrays
SPEC_VAL=         534.7675 / [nm] Mean wavelength
SPEC_BW =          312.949 / [nm] Bandpass width
TUTYP1  = 'Spectrum.Data.SpectralAxis.Value' /         
TUCD1   = 'em.wl;obs.atmos' / Air wavelength
TUNIT1  = 'angstrom'       /         
TCOMM1  = 'Computed from original WCS information' /         
TDMIN1  =          3782.93 /         
TDMAX1  =          6912.42 /         
TUTYP2  = 'Spectrum.Data.FluxAxis.Value' /         
TUCD2   = 'phot.flux.density;em.wl;stat.uncalib' /         
TUNIT2  = 'adu'            /         
TCOMM2  = 'Converted from 1-d pipeline spectrum (s1d_A)' /         
TDMIN2  = -96.20699999999999 /         
TDMAX2  =       2662.79858 /         
TUTYP3  = 'Spectrum.Data.FluxAxis.Accuracy.StatError' /         
TUCD3   = 'stat.error;phot.flux.density;em.wl;stat.uncalib' /         
TUNIT3  = 'adu'            /         
TCOMM3  = 'Error spectrum not available, filled with NaN' /         
CHECKSUM= 'ZFZEZ9YDZCYDZ9YD' / HDU checksum updated 2025-10-24T03:00:32
DATASUM = '2837289640'     / data unit checksum updated 2025-10-24T03:00:32
END
```

---

## 5. Summary of Best Practices
- Follow the **FITS keyword requirements** for wavelength calibration, flux calibration, and data quality.
- Use **Multi-Extension FITS (MEF)** for spectral data with multiple layers.
- Ensure metadata consistency and validation before submission.
- Utilize **ESO-provided validation tools** to check compliance.

---

> Note: For further details, refer to the official [ESO Phase 3 Documentation](https://www.eso.org/sci/observing/phase3/Phase3_doc.html) and the [Data Interface Control Document (DICD)](https://www.eso.org/sci/observing/phase3/DICD/DICD_latest.pdf).
