# IMAGE.FLUXMAP.md

## 1) What “Phase 3 compliant” means (quick context)

ESO Phase 3 is the process and standard for preparing, validating, and publishing **Science Data Products (SDPs)** in the ESO Science Archive. SDPs must be (a) calibrated in physical units, (b) traceable back to their provenance, and (c) richly described by FITS header keywords so the archive can index, validate, and serve them consistently. The **data category** you must write in the header for a flux map is:

```
PRODCATG = 'SCIENCE.IMAGE.FLUXMAP'
```

Flux maps are single 2‑D images with the science array in the **Primary HDU**. They represent **spectral flux density** per beam (e.g., `Jy/beam` or `mJy/beam`) and carry specific quality/astrometry parameters and APEX‑specific fields. Use **CD matrix WCS** (not CDELT) and include checksums.

---

## 2) Data definition (what the file must represent)

* **Science array**: 2‑D flux density in **Jansky per beam** (or decimal fractions) as declared by `BUNIT`.
* **Absolute flux scale**: Provide the overall **fractional calibration uncertainty** in `%` via `FLUXERR` (0–100).
* **Spatial resolution**: `SKY_RES` is the **FWHM effective beam** (arcsec); `SKY_RERR` = its RMS uncertainty (arcsec).
* **Astrometric uncertainty**: Use WCS error keywords `CSYERi` (systematic) and `CRDERi` (random).
* **Noise metric**: Use `BNOISE` for the **1σ point-source sensitivity** in **Jy** (for flux maps: equivalent to beam-smoothed RMS in Jy/beam). **Do not** use `ABMAGLIM` for flux maps.

---

## 3) Data format (how to store it)

* **Primary HDU contains the science image** (`NAXIS=2`). Flux maps are not MEF images.
* Use a **CDn_m matrix** for scale + rotation; **CDELTn not allowed** in Phase 3 ingestion. Set `CUNIT1/2='deg'`.
* **Filename & keywords**:

  * File suffix: `.fits` (uncompressed) or `.fits.fz` (fpack). Max header value length 68 characters; `LONGSTR`/`CONTINUE` not supported.
  * Include `CHECKSUM` & `DATASUM` in **all** HDUs.

---

## 4) Required associations (ancillary files)

* **Mandatory association**: the flux map **must** be associated to **RMS noise** and/or **SNR** map(s) using `ASSOCi`/`ASSONi`. *(Legend #33)*
* Recommended extras when applicable:

  * **EXPMAP** (exposure time per pixel, s) if exposure varies; `EXPTIME` (header) should be the median of EXPMAP values >0.
  * **GAINMAP**, **WEIGHTMAP/CONF**, **MASK** as needed.
* How to declare associations:

  * `ASSOCi = 'ANCILLARY.RMSMAP'` (or `SNRMAP`, etc.)
  * `ASSONi = 'filename.fits'`
  * `ASSOMi` (MD5) for **non-FITS** associated files. *(Legend #31)*

---

## 5) Header keyword tables

### 5.1 Primary HDU — **IMAGE.FLUXMAP**

| Keyword                                         |      Req. | Format / Typical value          | Notes                                                                      |
| ----------------------------------------------- | --------: | ------------------------------- | -------------------------------------------------------------------------- |
| `PRODCATG`                                      |    **MM** | `'SCIENCE.IMAGE.FLUXMAP'`       | Category marker.                                                           |
| `ASSOCi`                                        |  **MA31** | e.g. `'ANCILLARY.RMSMAP'`       | Use with `ASSONi`; required for non-FITS via `ASSOMi`. *(Legend #31)*      |
| `ASSONi`                                        |   **M33** | `'RMSMAP.fits'`/`'SNRMAP.fits'` | Flux maps **must** be associated to RMS and/or SNR map(s). *(Legend #33)*  |
| `ASSOMi`                                        |  **MA31** | MD5 string                      | For associated **non-FITS** files. *(Legend #31)*                          |
| `ORIGIN`                                        |    **MM** | Archive/institute string        | Data provider.                                                             |
| `TELESCOP`                                      |    **MM** | e.g. `APEX`                     | Telescope identifier.                                                      |
| `INSTRUME`                                      |    **MM** | e.g. `LABOCA`/`ARTEMIS`         | Instrument name.                                                           |
| `FILTER`/`FILTERi`                              |    **MM** | Instrument band ID              | Band used; instrument‑specific naming.                                     |
| `OBJECT`                                        |    **MM** | Target name                     | —                                                                          |
| `RA`, `DEC`                                     |    **MM** | degrees (J2000)                 | RA/DEC in footprint checks apply.                                          |
| `EQUINOX`                                       |  **MA37** | `2000.0`                        | Mandatory if `RADESYS='FK5'`; tolerated (2000.0) if `ICRS`. *(Legend #37)* |
| `RADESYS`                                       |    **MM** | `'ICRS'` or `'FK5'`             | Use with `EQUINOX` rule above.                                             |
| `TIMESYS`                                       |  **MA38** | e.g. `'UTC'`                    | Include if using a system **other than UTC**. *(Legend #38)*               |
| `EXPTIME`, `TEXPTIME`                           |    **MM** | seconds                         | For mosaics use EXPMAP median; `TEXPTIME` = total integration time.        |
| `MJD-OBS`, `MJD-END`                            |    **MM** | MJD (float)                     | `MJD-END ≥ MJD-OBS`.                                                       |
| `PROG_ID` / `PROGIDi`                           |    **ME** | programme IDs                   | Multi-programme rules apply.                                               |
| `OBIDi`                                         |    **ME** | OB identifiers                  | —                                                                          |
| `NCOMBINE`                                      |    **ME** | integer                         | Number of combined frames.                                                 |
| `OBSTECH`                                       |    **MM** | technique string                | E.g., mapping mode; consistent with instrument usage.                      |
| `FLUXCAL`                                       |    **MM** | `'ABSOLUTE'`                    | Imaging-type SDPs must be flux‑calibrated; **UNCALIBRATED** not allowed.   |
| `PROCSOFT`                                      |    **MM** | `name version`                  | Processing software.                                                       |
| `REFERENC`                                      |   **M40** | DOI/ADS bibcode or empty        | Empty allowed until publication. *(Legend #40)*                            |
| `PROVi` / `PROVXTN`                             |    **ME** | list / `T`                      | Set `PROVXTN=T` if you include a `PHASE3PROVENANCE` extension.             |
| `BUNIT`                                         |    **MM** | `'Jy/beam'` or `'mJy/beam'`     | Flux density per beam; declare the unit **here**.                          |
| `CRVALi`, `CRPIXi`, `CTYPEi`, `CUNITi`, `CDi_j` |    **MM** | WCS keywords                    | Use **CD** matrix; `CUNIT1/2='deg'`. No `CDELTn`.                          |
| `CSYERi` / `CRDERi`                             |    **RC** | degrees                         | WCS systematic / random astrometric errors.                                |
| `FLUXERR`                                       |    **MM** | `%` (0–100)                     | Overall flux scale error.                                                  |
| `WAVELMIN`, `WAVELMAX`                          |    **MM** | `nm`                            | Physical passband limits of the map; be consistent with instrument band.   |
| `RA_ERR`, `DEC_ERR`                             |    **NO** | —                               | Not allowed for flux maps.                                                 |
| `BNOISE`                                        |    **MM** | Jy (or Jy/beam eq.)             | 1σ point-source limit; for flux maps equals beam‑smoothed RMS.             |
| `MAPMODE`                                       |    **MM** | e.g. `'OTF'`, `'SPIRALRAS'`     | APEX mapping modes (comma‑separated).                                      |
| `FEBEi`                                         |    **MM** | e.g. `'LABOCA-ABBA'`            | APEX frontend/backend combination.                                         |
| `STOKES`                                        |    **SU** | `'/I/'` etc.                    | Contact Phase 3 before submitting polarimetry.                             |
| `ABMAGLIM`                                      |    **NO** | —                               | Not applicable; use `BNOISE`.                                              |
| `SKY_RES`                                       |    **MM** | arcsec                          | Effective beam FWHM of the map.                                            |
| `SKY_RERR`                                      |  **MA50** | arcsec                          | RMS error of `SKY_RES` *(Legend #50)*.                                     |
| `ARCFILE`, `CHECKSUM`, `DATASUM`                | **RE/MM** | strings                         | Checksums required; `ARCFILE` written/overwritten by archive.              |
| `ORIGFILE`, `P3ORIG`                            |    **RE** | strings                         | Reserved; do not set in extensions.                                        |
| `NOESODAT`                                      |  **MA59** | `T/F`                           | Use for products from **non‑ESO facilities**. *(Legend #59)*               |

> **Requirement code legend:** MM=Mandatory, ME=Mandatory (ESO), MA=Mandatory when applicable, RC=Recommended, NO=Not allowed, RE=Reserved, SU=Ask support. (See legend notes referenced above.)

---

### 5.2 Extension HDU(s) — minimal expectations

Flux maps **store science data in the Primary HDU**. If you add auxiliary HDUs (e.g., provenance table), keep them consistent:

| Keyword               |     Req. | Format / Typical value    | Notes                                                               |
| --------------------- | -------: | ------------------------- | ------------------------------------------------------------------- |
| `EXTNAME`             | **MA65** | e.g. `'PHASE3PROVENANCE'` | Specific values required in some cases (provenance, catalog links). |
| `ARCFILE`             |   **RE** | string                    | Archive will set/overwrite.                                         |
| `CHECKSUM`, `DATASUM` |   **MM** | strings                   | Include in **every** HDU.                                           |
| `ORIGFILE`, `P3ORIG`  |   **RE** | strings                   | Reserved by ESO archive.                                            |

If you include separate **ERROR/QUALITY images as extensions** (rare for flux maps), follow **HDUCLASn** and cross‑linking (`SCIDATA`, `ERRDATA`, `QUALDATA`).

---

## 6) WCS, units, and time — do‑this rules

* **Sky WCS**: 2‑axis RA/DEC with `CUNIT1='deg'`, `CUNIT2='deg'`, **CD matrix** only; set `CTYPEi` to a valid projection (e.g., `RA---TAN`, `DEC--TAN`).
* **Astrometry errors**: Provide `CSYERi` (systematic) and `CRDERi` (random) if known.
* **Flux units**: `BUNIT='Jy/beam'` (or `mJy/beam`). Note `BNOISE` is **always in Jy (or Jy/beam equivalent)** and relates to AB magnitude via the standard formula.
* **Times**: `MJD-OBS`, `MJD-END` float days; include `TIMESYS` **only if** not UTC. Keep `EXPTIME` consistent with the EXPMAP definition.

---

## 7) Validation expectations (what the Phase 3 checker looks for)

* **General**: header value lengths ≤ 68; **no** `LONGSTR`/`CONTINUE`; checksums present.
* **Coordinates**: `CUNIT1/2='deg'`; RA/DEC footprint sanity checks. **No `CDELTn`.**
* **Associations**: For flux maps, presence of **RMSMAP/SNRMAP** association is enforced.
* **Flux scale**: `FLUXCAL='ABSOLUTE'`; `FLUXERR` within 0–100.
* **Forbidden**: `ABMAGLIM` in flux maps. Use `BNOISE`.

---

## 8) Minimal header skeleton (illustrative)

```fits
SIMPLE  = T
BITPIX  = -32
NAXIS   = 2
PRODCATG= 'SCIENCE.IMAGE.FLUXMAP'
ORIGIN  = 'Your-Institute'
TELESCOP= 'APEX'
INSTRUME= 'LABOCA'
FILTER  = '870um'
OBJECT  = 'TargetName'
RADESYS = 'ICRS'
EQUINOX = 2000.0
CRVAL1  = ...
CRVAL2  = ...
CRPIX1  = ...
CRPIX2  = ...
CTYPE1  = 'RA---TAN'
CTYPE2  = 'DEC--TAN'
CUNIT1  = 'deg'
CUNIT2  = 'deg'
CD1_1   = ...
CD1_2   = ...
CD2_1   = ...
CD2_2   = ...
CSYER1  = ...
CRDER1  = ...
CSYER2  = ...
CRDER2  = ...
BUNIT   = 'mJy/beam'
FLUXCAL = 'ABSOLUTE'
FLUXERR = 12.5
SKY_RES = 19.5
SKY_RERR= 1.0
BNOISE  = 0.45
MAPMODE = 'OTF'
FEBE1   = 'LABOCA-ABBA'
EXPTIME = 3600.
TEXPTIME= 5400.
MJD-OBS = 60200.1234
MJD-END = 60200.1899
PROG_ID = '0101.A-0123'
OBID1   = '1234567'
NCOMBINE= 12
OBSTECH = 'IMAGE,OTF'
WAVELMIN= 350000.   / nm (example for submm band center-range)
WAVELMAX= 390000.
ASSOC1  = 'ANCILLARY.RMSMAP'
ASSON1  = 'my_rmsmap.fits'
REFERENC= ''         / or DOI/ADS id
CHECKSUM= '...'      / set by writer
DATASUM = '...'
END
```

> Replace `...` with real values; keep units as shown. Use instrument‑specific band limits for `WAVELMIN/MAX` in **nm**.

---

## 9) Submission checklist (copy/paste)

* [ ] `PRODCATG='SCIENCE.IMAGE.FLUXMAP'`
* [ ] Primary HDU holds the **science image** (no MEF image stack).
* [ ] WCS uses **CD** matrix; `CUNIT1/2='deg'`; no `CDELTn`.
* [ ] `BUNIT` is **Jy/beam** (or mJy/beam).
* [ ] `FLUXCAL='ABSOLUTE'`; `FLUXERR` is a % (0–100).
* [ ] `SKY_RES` (arcsec) and (if applicable) `SKY_RERR`.
* [ ] `BNOISE` present; **do not** include `ABMAGLIM`.
* [ ] `ASSOCi`/`ASSONi` link to **RMSMAP** and/or **SNRMAP**.
* [ ] APEX fields set if applicable: `MAPMODE`, `FEBEi`.
* [ ] RA/DEC, `RADESYS` (+ `EQUINOX` rule), `EXPTIME`/`TEXPTIME`, `MJD-OBS/END`.
* [ ] `PROG_ID`/`OBIDi`/`NCOMBINE`, `OBSTECH`, `PROCSOFT`, `REFERENC`.
* [ ] `CHECKSUM`/`DATASUM` in **all** HDUs.
* [ ] Reserved keys (`ARCFILE`, `ORIGFILE`, `P3ORIG`) **not manually set** (archive will manage them).

---

## 10) Ambiguities & caveats

* **`WAVELMIN/WAVELMAX` in flux maps**: Not explicitly defined like for spectra/cubes; set to the **physical passband limits** (nm) of the instrument/filter used for the map so archive services can index wavelength coverage. Keep consistent with instrument band definition and header units.
* **Polarimetry (`STOKES`)**: Contact Phase 3 operations **before** submitting polarized flux maps; additional conventions apply.

---

## 11) References (quick pointers)

* **SDP standard**: filename rules; categories & associations; HDUCLASn & checksums; `BNOISE`; APEX (`MAPMODE`, `FEBEi`); keyword matrices & appendices.
* **DIC**: WCS and file‑structure guidance; prefer **CD matrix** in SDPs.
* **Legend & checks**: Keyword requirement codes and Phase 3 checker expectations.