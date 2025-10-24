# **SPECTRUM** – Phase 3 Science Data Products (SDP) Submission Guide

## 1. Introduction

**Purpose.** This guide condenses the ESO SDP standard and DICD into clear, practical steps for preparing and submitting **SCIENCE.SPECTRUM** products (1D spectra) to Phase 3.

**What counts as SPECTRUM.** A **single-target 1D spectrum** delivered as a **single-row FITS BINTABLE** with all required metadata and consistent arrays.

**Key principles.**

* Products are fully calibrated, documented, and traceable.
* FITS headers capture observation context, calibration status, and VO interoperability.

### 1.1 Overview of how to submit data to Phase 3

1. Prepare FITS files per this guide (structure + keywords).
2. Provide a concise release description (methods, calibrations, data quality).
3. Validate locally (headers, VO fields, checksums) and then with Phase 3 tools.
4. Deliver via the Phase 3 process (batch with unique filenames; include ancillary files if needed).

---

## 2. Spectral Data Submission Requirements

### 2.1 General Data Format (SPECTRUM)

* **Primary HDU:** header only (**NAXIS=0**).
* **One extension only:** **BINTABLE** with **NAXIS=2**, **NAXIS2=1**, **PCOUNT=0**, **GCOUNT=1**.
* **Arrays & fields:**

  * **Mandatory columns (TTYPEi):** spectral coordinate (**WAVE** or **FREQ** or **ENER**), **FLUX**, **ERR** (in that order).
  * **Optional columns:** **QUAL**, **BGFLUX**, **CONTINUUM**, **EXPOSURE** (and others if justified).
  * **All arrays same length:** **NELEM** (extension header) and consistent **TFORMi**. No variable-length arrays.
* **Units:** Do **not** use `BUNIT` for 1D spectra; specify per-column **`TUNITi`**. If dimensionless, set `TUNITi=''`.
* **EXTNAME:** every extension name must be unique within the file.

### 2.2 Wavelength Calibration

* Spectral axis may be **wavelength**, **frequency**, or **energy**.
* Use **`WAVE`** with TUCD `em.wl` (vacuum) or `em.wl;obs.atmos` (air). State one consistently.
* Provide spectral coverage with **`WAVELMIN`**/**`WAVELMAX`**.
* For wavelength arrays, values **strictly increase**; define **`TDMIN1`/`TDMAX1`** for start/stop.

### 2.3 Flux Calibration and Photometry Requirements

* Set **`FLUXCAL`** to describe calibration status. If **absolute** flux calibrated:

  * Provide **`FLUXERR`** (global absolute-flux uncertainty) and ensure **`TUNIT(FLUX)`** is a physical unit.
* If serving **continuum-normalised** flux as default:

  * **`CONTNORM = T`** in primary; main flux column is normalised (set `TUNIT=''`, and use UCD with `arith.ratio`).
* If providing multiple flux versions (e.g., recipe A/B), use **VOCLASS='SPECTRUM v2.0'** and mark the **main** flux/error with `meta.main` in **`TUCDi`**.

### 2.4 Temporal Information

* Provide **`MJD-OBS`** (start), **`MJD-END`** (end), and optionally **`TMID`** (midpoint) in the extension.
* Include **`EXPTIME`** and **`TEXPTIME`** in primary; if you supply an **`EXPOSURE`** array column, set its `TUCD` to `time.duration;obs.exposure`.

---

## 3. Required and Recommended FITS Header Keywords

> The tables below list the **minimum** metadata for SPECTRUM products. “Mandatory (ESO)” denotes keywords required in ESO context; keep them if available in your provenance. “Conditional” indicates rules in the Notes beneath.

### 3.1 Primary HDU Keywords

| Keyword                | Description                                                     | Requirement               |
| ---------------------- | --------------------------------------------------------------- | ------------------------- |
| `PRODCATG`             | Data product category; set to `SCIENCE.SPECTRUM`.               | Mandatory                 |
| `ORIGIN`               | Data origin (e.g., facility/site).                              | Mandatory                 |
| `TELESCOP`             | Telescope identifier.                                           | Mandatory                 |
| `INSTRUME`             | Instrument identifier.                                          | Mandatory                 |
| `OBJECT`               | Target name (well-established name for moving objects allowed). | Mandatory                 |
| `RA`, `DEC`            | Target coordinates (ICRS/FK5).                                  | Mandatory                 |
| `EQUINOX`              | Epoch (required if `RADESYS='FK5'`, value 2000.0).              | Conditional               |
| `RADESYS`              | Coordinate reference frame (`ICRS`, `FK5`).                     | Mandatory                 |
| `TIMESYS`              | Time system (`UTC` or `TAI` only). If not UTC, specify.         | Conditional               |
| `MJD-OBS`, `MJD-END`   | Observation start/end (MJD).                                    | Mandatory                 |
| `EXPTIME`, `TEXPTIME`  | Exposure time(s) [s].                                           | Mandatory                 |
| `PROG_ID`              | Programme ID.                                                   | Mandatory (ESO)           |
| `OBIDi`                | Observation Block ID(s).                                        | Mandatory (ESO)           |
| `NCOMBINE`             | Number of combined exposures.                                   | Mandatory (ESO)           |
| `OBSTECH`              | Observation technique (e.g., SPECTRUM, ECHELLE, MOS).           | Mandatory                 |
| `FLUXCAL`              | Flux calibration status.                                        | Mandatory                 |
| `PROCSOFT`             | Processing software (name/version).                             | Mandatory                 |
| `REFERENC`             | Bibliographic reference (may be empty if no paper yet).         | Mandatory (note)          |
| `PROVi` / `PROVXTN`    | Processing provenance keywords / extension.                     | Mandatory (ESO)           |
| `SPECSYS`              | Spectral reference frame (e.g., BARYCENT, HELIOCEN…).           | Mandatory                 |
| `CONTNORM`             | `T` if continuum-normalised spectrum is served as default.      | Mandatory when applicable |
| `TOT_FLUX`             | Total integrated flux (if relevant to product).                 | Mandatory when applicable |
| `FLUXERR`              | Absolute flux calibration uncertainty (rules in Notes).         | Conditional               |
| `WAVELMIN`, `WAVELMAX` | Spectral coverage (consistent with air/vacuum choice).          | Mandatory                 |
| `SPEC_BIN`             | Spectral bin width or sampling.                                 | Mandatory                 |
| `SPEC_ERR`, `SPEC_SYE` | Random/systematic spectral calibration errors.                  | Recommended               |
| `RA_ERR`, `DEC_ERR`    | Positional uncertainties (1D spectra).                          | Recommended               |
| `EXT_OBJ`              | `T` for extended sources, else `F`.                             | Mandatory when applicable |
| `CHECKSUM`, `DATASUM`  | FITS checksums (primary and extension HDUs).                    | Mandatory                 |

**Notes (Primary).**

* `EQUINOX` required **only** with `RADESYS='FK5'` (value 2000.0).
* If **absolute** flux calibrated: set `FLUXERR` to the global flux-scale uncertainty; if uncalibrated or not determinable, use sentinel values as defined by standard practice.
* `BUNIT` is **not used** in SPECTRUM primary HDU; use `TUNITi` per column instead.

### 3.2 Extension HDU Keywords

| Keyword            | Description                                                 | Requirement                 |
| ------------------ | ----------------------------------------------------------- | --------------------------- |
| `NELEM`            | Pixel length of data arrays (all columns equal length).     | Mandatory                   |
| `VOCLASS`          | VO data model (`SPECTRUM v1.0` or `SPECTRUM v2.0`).         | Mandatory                   |
| `VOPUB`            | Publisher (`ESO/SAF`).                                      | Mandatory                   |
| `TITLE`            | Short, human-readable dataset title (unique in collection). | Mandatory                   |
| `APERTURE`         | Slit width or fibre diameter [deg].                         | Mandatory                   |
| `TELAPSE`          | Total elapsed time `MJD-END − MJD-OBS` [s].                 | Mandatory                   |
| `TMID`             | Exposure midpoint `(MJD-OBS + MJD-END)/2`.                  | Mandatory                   |
| `SPEC_VAL`         | Characteristic spectral coordinate [(WMAX+WMIN)/2] [nm].    | Mandatory                   |
| `SPEC_BW`          | Spectral width `(WMAX − WMIN)` [nm].                        | Mandatory                   |
| `TTYPEi`           | Column labels (e.g., `WAVE`, `FLUX`, `ERR`, `QUAL`…).       | Mandatory                   |
| `TFORMi`           | Column formats (FITS BINTABLE).                             | Mandatory                   |
| `TUNITi`           | Column units (empty string if dimensionless).               | Mandatory                   |
| `TUTYPi`           | VO Utype (use IVOA Spectrum DM; set `''` if undefined).     | Mandatory                   |
| `TUCDi`            | UCDs (mark main flux/error with `meta.main` when used).     | Mandatory                   |
| `TDMIN1`, `TDMAX1` | Start/stop spectral coordinate for column 1.                | Mandatory for spectral axis |
| `EXTNAME`          | Unique per file; descriptive name for the table HDU.        | Mandatory when applicable   |

**Mandatory columns (order)**

1. `WAVE`/**`FREQ`**/**`ENER`**  → spectral coordinate
2. `FLUX` (or `FLUX_*`) → main data spectrum
3. `ERR` (or `ERR_*`) → error spectrum

**Common optional columns**

* `QUAL` (integer data quality mask), `BGFLUX` (sky spectrum, at most one), `CONTINUUM` (continuum model), `EXPOSURE` (per-pixel exposure time).

**VO mapping & multiple fluxes**

* Single flux: `VOCLASS='SPECTRUM v1.0'`, use IVOA utypes.
* Multiple fluxes (e.g., `FLUX_A`,`FLUX_B`): `VOCLASS='SPECTRUM v2.0'`; tag the default flux/error with `meta.main` in `TUCD` and use consistent `TUTYP` prefixes (e.g., `spec:` vs `eso:`).

### 3.3 Additional Considerations

* **EXTNAME uniqueness** across all HDUs in the file.
* **Telluric-corrected flux** columns can be named `FLUX_TELLURIC`.
* For **air vs. vacuum wavelengths**, ensure `WAVELMIN/MAX` match the spectral-axis choice and `TUCD1` reflects air (`em.wl;obs.atmos`) or vacuum (`em.wl`).

---

## 4. File Structure and Validation

**Structure checklist**

* Primary HDU: header only (`NAXIS=0`).
* One **BINTABLE** extension with `NAXIS2=1`, `PCOUNT=0`, `GCOUNT=1`.
* Arrays consistent length; **no** variable-length arrays.
* `CHECKSUM`/`DATASUM` present in primary and extension.

**Validation checklist**

* Header passes FITS validation; values/units agree with VO and Phase 3 rules.
* `RA/DEC`, `MJD-OBS/END`, `SPECSYS`, `WAVELMIN/MAX`, `FLUXCAL` present and consistent.
* `VOCLASS`, `TTYPE/TFORM/TUNIT/TUCD/TUTYP` defined; `TDMIN1/TDMAX1` for spectral axis.
* For continuum-normalised products: `CONTNORM=T`, `TUNIT(FLUX)=''`, and UCDs include `arith.ratio`.
* EXTNAMEs unique; if provenance or associated data are stored in auxiliary HDUs/files, cross-reference cleanly.

---

## 5. Summary of Best Practices

* Keep the file minimal: **Primary header only**, a **single** BINTABLE extension.
* Use **standard column names** (`WAVE/FREQ/ENER`, `FLUX`, `ERR`, optional `QUAL`, `BGFLUX`, `CONTINUUM`, `EXPOSURE`).
* **Be explicit** about **air vs. vacuum** wavelengths and **spectral frame** (`SPECSYS`).
* **Mark the main flux/error** with `meta.main` in `TUCD` when multiple flux versions exist; set `VOCLASS` accordingly.
* Always provide **`WAVELMIN/MAX`** and **`SPEC_VAL/BW`** consistent with the spectral axis.
* **No `BUNIT`** in SPECTRUM; use `TUNITi` per column.
* Include **processing provenance** and **software version**; keep **checksums**.
* Validate with Phase 3 tools before delivery; ensure **filenames are unique** within the batch.

---

### Appendix – Quick Column/UCD/Unit hints

* `WAVE` → UCD `em.wl` (vacuum) or `em.wl;obs.atmos` (air); `TUNIT='nm'` (or suitable).
* `FLUX` (absolute) → UCD `phot.flux.density;em.wl` (add `;meta.main` if default); typical `TUNIT='erg s-1 cm-2 Å-1'` (or SI).
* `FLUX_NORM` (continuum-normalised) → include `arith.ratio` in UCD; `TUNIT=''`.
* `ERR` → `stat.error;phot.flux.density;em.wl` (add `;meta.main` if default).
* `EXPOSURE` → `time.duration;obs.exposure`; `TUNIT='s'`.
