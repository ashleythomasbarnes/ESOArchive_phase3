#!/usr/bin/env python
"""
APEX Phase 3 Tool: Creation of FITS spectral cubes compliant with ESO Phase 3 requirements.
Contact: pvenegas@eso.org, tstanke@mpe.mpg.de
Creation: 2020-07-13; Updated: 2021-08-13
"""

import os
import sys
import math
import copy
import hashlib
import datetime
import pathlib
import re

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from spectral_cube import SpectralCube
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs

# Colors for printing messages
RED = "\033[31m"       # For errors
GREEN = "\033[32;3m"   # For outputs
BLUE = "\033[34;1m"    # Request inputs
MAGENTA = "\033[35m"   # Header
CYAN = "\033[36m"      # Main words
ENDCOLOR = "\033[0m"

##########################################
# CONFIGURATION: Set input parameters here.
##########################################
# Set the names of the input FITS files (must be present in the current directory)
SPEC_IN = "G030_13CO21_Tmb_DR1.fits"           # Science cube file
RMS_IN = "G030_13CO21_Tmb_DR1.fits"            # Error (rms) cube file
FILENAME_OUT = SPEC_IN.split('_Tmb_DR1.fits')[0]+'_13CO21_P3.fits'
FILENAME_WHITE_OUT = SPEC_IN.split('_Tmb_DR1.fits')[0]+'_13CO21_P3_whitelight.fits'

# LMV file association (if applicable)
ASSOCIATE_LMV = False                       # Set True if you want to associate a .lmv file
LMV_FILE = "your_cube.lmv"                  # Specify the .lmv file (if ASSOCIATE_LMV is True)

# TAP file pause: set to True if you wish to pause for manual TAP file editing
TAP_PAUSE = False

# Jy/K conversion factor handling:
# If DIFFERENT_JY_FACTOR is True, the default factor from the dictionary will be overridden.
DIFFERENT_JY_FACTOR = False
NEW_JY_FACTOR = 53.0                        # New Jy/K factor (if overriding)
NEW_JY_UNCERT = 8.0                         # New Jy/K uncertainty

# Beam efficiency correction:
# If APPLY_BEAM_EFFICIENCY is True, the beam efficiency correction will be applied.
APPLY_BEAM_EFFICIENCY = False
BEAM_EFFICIENCY = 1.0                        # Value between 0 and 1

# Flux conversion:
# Set CONVERT_TO_JY to True if you want to convert from K to Jy.
CONVERT_TO_JY = False

# Programme Identification:
# If CONFIRM_PROG_CODE is False, the new programme code NEW_PROG_CODE will be used.
CONFIRM_PROG_CODE = True
NEW_PROG_CODE = ""                          # Only used if CONFIRM_PROG_CODE is False

# Bibliographic reference (Bibcode); leave empty if none.
BIBLIO_REF = "2021MNRAS.500.3064S"

# Beam size override:
# Set NEW_BEAM_SIZE to a number (in arcsec) to override the effective beam size;
# otherwise leave as None to use the value computed from the header.
NEW_BEAM_SIZE = None

##########################################
# Dictionaries for Jy/K conversion factors
##########################################
dicc_Gen = {
    "HET230": [2016, 39, 6],
    "PI230": [2016, 44, 6],
    "HET345": [2016, 53, 8],
    "NOVA660": [2016, 137, 22],
    "GARD180": [2016, 39, 6],
    "PI230": [2017, 45, 7],
    "HET230": [2017, 40, 6],
    "HET345": [2017, 51, 8],
    "NOVA660": [2017, 110, 18],
    "GARD180": [2017, 40, 6],
    "HET460": [2017, 72, 11],
    "PI230": [2018, 46, 4],
    "SEPIA660": [2018, 63, 5],
    "GARD180": [2018, 40, 6],
    "SEPIA180": [2019, 35, 3],
    "SEPIA660": [2019, 68, 6],
}

dicc = {
    "HET230": [2016, 40, 6],
    "PI230": [2017, 45, 7],
    "HET345": [2016, 53, 8],
    "SUPERCAM": [2014, 53, 8],  # Use same value as for FLASH345 and HET345 also for SUPERCAM
    "FLASH345": [2018, 53, 4],
    "HET460": [2017, 72, 11],
    "GARD180": [2017, 40, 6],
    "NOVA660": [2017, 110, 18],
    "SEPIA180": [2019, 35, 3],
    "SEPIA660": [2019, 68, 6],
}

##########################################
# Function definitions
##########################################
def tapQuery_Het(source):
    """
    Query the ESO TAP service.
    """
    from pyvo.dal import tap
    from astropy.coordinates import SkyCoord
    from astropy.units import Quantity
    from tabulate import tabulate

    ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
    tapobs = tap.TAPService(ESO_TAP_OBS)

    long = np.int32(source.split('G')[1])
    source_left = f'G{long-1}'
    source_right = f'G{long+1}'

    print("\nQuerying the ESO TAP service at %s" % ESO_TAP_OBS)
    query = (
        "SELECT dp_id, exposure, prog_id, object, dp_tech, instrument, ra, dec \n"
        "FROM dbo.raw \n"
        "WHERE dp_id LIKE 'APEXHET.%' \n"
        "AND (prog_id LIKE '092.F-9315%' OR prog_id LIKE '193.C-0584%') \n"
        f"AND (object LIKE '{source}%' OR object LIKE '{source_left}%' OR object LIKE '{source_right}%') \n"
        "AND dp_cat = 'SCIENCE'"
    )
    print("\nQuery:\n" + query)
    res = tapobs.search(query=query, maxrec=1000)
    print("\n", res.to_table(), "\n")
    print("\nA total of " + str(len(res.to_table())) + " records were found matching the provided criteria.")
    table = res.to_table()
    filename = f"{source}.tap"
    print("\nResults has been written to: " + filename)
    with open(filename, "w") as f:
        print(tabulate(table), file=f)
    return table


def get_frontBack():
    """
    Get instrument-backend information (currently hardcoded).
    """
    # In this version, FEBE is set via the TAP query section below.
    return "SUPERCAM-SCBE"


def freqArray(prihdr, specNr):
    """
    Compute frequency array from the WCS in the header.
    """
    v_light = 299792458.0  # Speed of light (m/s)
    cdelv = prihdr["CDELT3"]
    crpix1 = prihdr["CRPIX3"]
    crvel = prihdr["CRVAL3"]
    refreq = prihdr["RESTFREQ"]
    cdelf = -1.0 * refreq * cdelv / v_light
    crfreq = refreq * (1 - crvel / v_light)
    freq0 = crfreq - (crpix1 - 1) * cdelf
    freqs = []
    print("Frecuency axis in construction")
    for i in range(specNr):
        f1 = (freq0 + i * cdelf) / 1.0e9  # Convert to GHz
        freqs.append(f1)
    return freqs, refreq, cdelf

def check_if_fk5(hdu):
    """ Check if the WCS is in FK5 """
    if ('RA' in hdu.header['CTYPE1']) & ('DEC' in hdu.header['CTYPE2']):
        return True
    else:
        return False

def convert_fk5(hdu):
    """ Convert the WCS to FK5 """

    print('Converting to FK5...')

    hdr = hdu.header
    hdu_2d = hdu.copy()
    hdu_2d.data = hdu_2d.data[0, :, :]
    hdu_2d.header['NAXIS'] = 2
    hdu_2d.header['WCSAXES'] = 2
    del hdu_2d.header['*3*']

    wcs_out, shape_out = find_optimal_celestial_wcs(hdu_2d, frame='fk5', auto_rotate=True, projection='GLS')
    # wcs_out, shape_out = find_optimal_celestial_wcs(hdu_2d, frame='fk5')

    hdr_out_ = wcs_out.to_header()

    hdr_out = hdr.copy()
    for key in ['CRPIX1', 'CRPIX2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1', 'CDELT2', 
                'CUNIT1', 'CUNIT2', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'LONPOLE', 'LATPOLE', 
                'RADESYS', 'EQUINOX']:
        try: 
            hdr_out[key] = hdr_out_[key]
        except: 
            print('Key not found:', key)
            continue
    del hdr_out['*PV*']

    array_out = reproject_interp(hdu, hdr_out, return_footprint=False)

    hdu_out = fits.PrimaryHDU(array_out, header=hdr_out)

    cube = SpectralCube.read(hdu_out)
    cube.write('tmp.fits', overwrite=True)
    hdu = fits.open('tmp.fits')
    os.system('rm tmp.fits')

    return hdu

##########################################
# Main script execution
##########################################
def main():
    # --- Print header ---
    print("##########################################")
    print(MAGENTA + "#APEX-Phase 3 Tool" + ENDCOLOR)
    print(MAGENTA + "#Creation of FITS spectral cubes" + ENDCOLOR)
    print(MAGENTA + "#compliant with ESO Phase 3 requirements" + ENDCOLOR)
    print(MAGENTA + "#ESO-APEX 12m telescope" + ENDCOLOR)
    print(MAGENTA + "#contact: pvenegas@eso.org, tstanke@mpe.mpg.de" + ENDCOLOR)
    print("##########################################\n")

    # --- Important information ---
    print(RED + "Important information:" + ENDCOLOR, "\n",
          "- All the RAW science data files (one per scan number) used to produce this data product",
          "should be present in the current folder (in ARCFILE format with extension .fits). Do not",
          "include calibration data.", "\n",
          "- RAW science files can be downloaded from http://archive.eso.org/wdb/wdb/eso/apex/form.",
          'Run the query to find your scan(s) and then make sure you select the "Mark RAW" column',
          "before requesting your data.\n")

    # --- Use configured input files ---
    print("\n" + CYAN + "Your current FITS files" + ENDCOLOR)
    os.system("ls *.fits")
    print("")

    listOfFiles = os.listdir(".")
    # Use configured science file
    specIn = SPEC_IN
    print("Using science FITS file: " + specIn)
    fits.info(specIn)
    specfile = fits.open(specIn)

    # The following block is commented out because the GAL conversion is accepted in the archive... 
    # # --- Check if the WCS is in FK5 ---
    # if ~check_if_fk5(specfile[0]):
    #     specfile = convert_fk5(specfile[0])

    # --- Open FLUX FITS file and extract key header/data info ---
    try:
        scidat = specfile[0].data
        prihdr = specfile[0].header
        orighdr = copy.copy(prihdr)
        print(prihdr[0])
        specNr = scidat.shape[0]
        b_unit = str(prihdr["BUNIT"])
        source = prihdr["OBJECT"]

        # Ensure RA/DEC keywords exist
        if ("RA" not in prihdr) and ("RA" in prihdr["CTYPE1"]):
            prihdr["RA"] = prihdr["CRVAL1"]
        if ("DEC" not in prihdr) and ("DEC" in prihdr["CTYPE2"]):
            prihdr["DEC"] = prihdr["CRVAL2"]
        else: 
            # Convert Galactic CRVAL1 and CRVAL2 to RA/DEC 
            print('Converting Galactic to RA/DEC')
            c = SkyCoord(l=orighdr['CRVAL1'], b=orighdr['CRVAL2'], frame='galactic', unit=(u.deg, u.deg))
            prihdr["RA"] = c.icrs.ra.deg
            prihdr["DEC"] = c.icrs.dec.deg

        ra = prihdr["RA"]
        dec = prihdr["DEC"]

        print("\nObject name : ", source, "; Flux unit: ", b_unit, "\n")
    except Exception:
        print(RED + "Wrong FLUX FITS file, try it once again." + ENDCOLOR)
        sys.exit()

    specfile.close()

    # --- Open error (rms) cube ---
    rmsIn = RMS_IN
    print("\nUsing error cube FITS file: " + rmsIn)
    if rmsIn not in listOfFiles:
        print(RED + f"file '{rmsIn}' containing error cube does not exist. Please provide a valid file." + ENDCOLOR)
        sys.exit()

    fits.info(rmsIn)
    rmsfile = fits.open(rmsIn)

    # # --- Check if the WCS is in FK5 ---
    # if ~check_if_fk5(rmsfile[0]):
    #     rmsfile = convert_fk5(rmsfile[0])

    # --- Open RMS FITS file and extract key header/data info ---
    print(CYAN + "... checking rms cube ..." + ENDCOLOR)
    try:
        rmsdat = rmsfile[0].data
        rmshdr = rmsfile[0].header
        b_unit_rms = str(rmshdr["BUNIT"])
        source_rms = rmshdr["OBJECT"]

        if ("RA" not in prihdr) and ("RA" in prihdr["CTYPE1"]):
            prihdr["RA"] = prihdr["CRVAL1"]
        if ("DEC" not in prihdr) and ("DEC" in prihdr["CTYPE2"]):
            prihdr["DEC"] = prihdr["CRVAL2"]
        else: 
            # Convert Galactic CRVAL1 and CRVAL2 to RA/DEC 
            print('Converting Galactic to RA/DEC')
            c = SkyCoord(l=orighdr['CRVAL1'], b=orighdr['CRVAL2'], frame='galactic', unit=(u.deg, u.deg))
            prihdr["RA"] = c.icrs.ra.deg
            prihdr["DEC"] = c.icrs.dec.deg

        ra = prihdr["RA"]
        dec = prihdr["DEC"]

        print("\nObject name : ", source_rms, "; Flux unit: ", b_unit_rms, "\n")
    except Exception:
        print(RED + "Wrong rms FITS file, try it once again." + ENDCOLOR)
        sys.exit()

    if rmshdr["NAXIS"] != 3:
        print(RED + "rms file seems not to be a cube... Please give a different file!" + ENDCOLOR)
        sys.exit()
    else:
        print(GREEN + "      ... rms file is a cube ..." + ENDCOLOR)

    if (rmshdr["NAXIS1"] != prihdr["NAXIS1"] or rmshdr["NAXIS2"] != prihdr["NAXIS2"] or rmshdr["NAXIS3"] != prihdr["NAXIS3"]):
        print(RED + "rms and science cube dimensions disagree!" + ENDCOLOR)
        sys.exit()
    else:
        print(GREEN + "      ... rms and science cube dimensions agree ..." + ENDCOLOR)

    if (rmshdr["CRVAL1"] != prihdr["CRVAL1"] or rmshdr["CRPIX1"] != prihdr["CRPIX1"] or rmshdr["CDELT1"] != prihdr["CDELT1"] or
        rmshdr["CRVAL2"] != prihdr["CRVAL2"] or rmshdr["CRPIX2"] != prihdr["CRPIX2"] or rmshdr["CDELT2"] != prihdr["CDELT2"]):
        print(RED + "rms and science cube X/Y coordinate definitions disagree!" + ENDCOLOR)
        print(RED + "Please check the input rms file." + ENDCOLOR)
        sys.exit()
    else:
        print(GREEN + "      ... rms and science cube X/Y coordinate definitions agree ..." + ENDCOLOR)

    if (rmshdr["CRVAL3"] != prihdr["CRVAL3"] or rmshdr["CRPIX3"] != prihdr["CRPIX3"] or rmshdr["CDELT3"] != prihdr["CDELT3"]):
        print(RED + "rms and science cube spectral axis definitions disagree!" + ENDCOLOR)
        sys.exit()
    else:
        print(GREEN + "      ... rms and science cube spectral axis definitions agree ..." + ENDCOLOR)

    rmsfile.close()

    # --- LMV file association ---
    if ASSOCIATE_LMV:
        attlmv = "y"
        lmvIn = LMV_FILE
    else:
        attlmv = "n"
        lmvIn = ""

    if attlmv == "y":
        if (os.path.exists(lmvIn)) and (lmvIn.endswith("lmv")):
            print(CYAN + f"file {lmvIn} exists; a copy will be associated to the Phase 3 cube" + ENDCOLOR)
        else:
            print(RED + f"file '{lmvIn}' does not exist or has wrong extension." + ENDCOLOR)
            sys.exit()
    else:
        print(GREEN + "No .lmv file will be associated, continuing" + ENDCOLOR)

    # --- TAP query ---
    queryfile = tapQuery_Het(source)
    print(CYAN + "If you need to modify the TAP file, the time is now!" + ENDCOLOR)
    if TAP_PAUSE:
        input("Please press the Enter key to resume\n")
    else:
        print("Skipping TAP file pause.\n")

    # Read the TAP file
    fileq = f"{source}.tap"
    with open(fileq, "r") as queryfile:
        qfile = [line.split() for line in queryfile][1:-1]
    print(BLUE + "After pause, your TAP file length is:" + ENDCOLOR)
    print("Length = " + str(len(qfile)) + " rows\n")
    fileList = []
    extCont = []
    proList = []
    sourceList_tmp = []
    dprtech_tmp = []
    instrum = []
    posra = []
    posdec = []
    for row in qfile:
        fileList.append(row[0])
        extCont.append(float(row[1]))
        proList.append(row[2])
        sourceList_tmp.append(row[3])
        dprtech_tmp.append(row[4])
        instrum.append(row[5])
        posra.append(row[6])
        posdec.append(row[7])
    fileList.sort()
    dprtech = list(set(dprtech_tmp))
    obs_tech = dprtech[0]
    instrument = instrum[0]
    print(CYAN + "These files will be associated with the cube:\n " + "\n".join(fileList) + ENDCOLOR)
    NrFilesBol = len(fileList)
    if NrFilesBol >= 1:
        print("\n" + GREEN + "Total <APEXHET.*.fits> files >> " + ENDCOLOR, NrFilesBol)
        print(*fileList, sep="\t")
        print("")
    else:
        print(RED + "No APEXHET files found in the current directory." + ENDCOLOR)
        sys.exit()

    proList = list(set(proList))
    if len(proList) == 1:
        prog_code = proList[0]
    elif len(proList) > 1:
        prog_code = "MULTI"
    print(prog_code)

    ext_val = sum(extCont)
    if ext_val < 5:
        print(RED + "Integration time inconsistency. Check your data" + ENDCOLOR)
        sys.exit()

    listTime = []
    for tmp_t in fileList:
        fileTime = tmp_t[8:18] + " " + tmp_t[19:31]
        listTime.append(fileTime)
        print(tmp_t, fileTime)

    # --- FEBE and Jy/K conversion factor ---
    febe = "SUPERCAM-SCBE"  # Hardcoded in this pipeline
    ind = febe.index("-")
    front = febe[:ind]
    fac_yr = dicc[front][0]
    factor = dicc[front][1]
    fac_un = dicc[front][2]
    print(obs_tech, febe, front, factor)
    print(CYAN + "APEX spectral line cubes may be ingested in the archive either on a \nbrightness temperature scale (K) or flux density (Jy) scale." + ENDCOLOR)
    print("")
    print(CYAN + "- To obtain the conversion factor (Jy/K), please refer to the factor that best fits the time range of your observations." + ENDCOLOR)
    print("http://www.apex-telescope.org/telescope/efficiency/index.php")
    print("")
    print(CYAN + f"Currently used Jy/K factor for {front}: {factor}+/-{fac_un}  (Year: {fac_yr})" + ENDCOLOR)
    print("")
    if DIFFERENT_JY_FACTOR:
        factor = NEW_JY_FACTOR
        fac_un = NEW_JY_UNCERT
        print(GREEN + "New Jy/K factor set to: " + ENDCOLOR, factor, "+/-", fac_un)
    else:
        print(GREEN + "Using default Jy/K factor" + ENDCOLOR)

    # --- Flux conversion ---
    try:
        if b_unit[0] == "K":
            print(CYAN + "- Your data product is currently on a brightness temperature scale [K]." + ENDCOLOR)
            # Assume no additional beam efficiency correction unless APPLY_BEAM_EFFICIENCY is True.
            if APPLY_BEAM_EFFICIENCY:
                beam_eff = BEAM_EFFICIENCY
                if 0.0 < beam_eff <= 1.0:
                    print(CYAN + f"... applying beam efficiency correction: {beam_eff}" + ENDCOLOR)
                    scidat = scidat / beam_eff
                    rmsdat = rmsdat / beam_eff
                else:
                    print(RED + "Beam efficiency must be between 0 and 1." + ENDCOLOR)
                    sys.exit()
            else:
                print(CYAN + "... no additional beam efficiency correction applied" + ENDCOLOR)
            if CONVERT_TO_JY:
                flux = scidat * factor
                flux_rms = rmsdat * factor
                b_unit = "Jy"
                b_unit_rms = "Jy"
            else:
                flux = scidat
                flux_rms = rmsdat
        elif b_unit == "mJy":
            print(CYAN + "- Your data product is in Flux density [mJy]" + ENDCOLOR)
            if CONVERT_TO_JY:
                flux = scidat / (1000.0 * factor)
                flux_rms = rmsdat / (1000.0 * factor)
                b_unit = "K"
                b_unit_rms = "K"
            else:
                flux = scidat
                flux_rms = rmsdat
        elif b_unit == "Jy":
            print(CYAN + "- Your data product is in Flux density [Jy]" + ENDCOLOR)
            if CONVERT_TO_JY:
                flux = scidat / factor
                flux_rms = rmsdat / factor
                b_unit = "K"
                b_unit_rms = "K"
            else:
                flux = scidat
                flux_rms = rmsdat
        else:
            print(RED + "Data not in recognized units (K, mJy, or Jy)." + ENDCOLOR)
            sys.exit()
    except Exception:
        print(RED + "Error in flux conversion. Please check intensity units." + ENDCOLOR)
        sys.exit()

    rms_med = np.nanmedian(flux_rms)
    rms_med_Jy = rms_med * factor if b_unit_rms[0] == "K" else rms_med
    flux_white = np.nanmean(flux, axis=0)
    freqs, refreq, cdel = freqArray(prihdr, specNr)

    # --- Create HDUs ---
    hdu_fluxcube = fits.ImageHDU(data=flux, header=None, name="DATA_EXT")
    hdu_rmscube = fits.ImageHDU(data=flux_rms, header=None, name="STAT_EXT")
    hdu_white = fits.PrimaryHDU(flux_white)

    # --- Observation timing ---
    time_1 = listTime[0]
    time_2 = listTime[-1]
    tmp_t1 = Time(time_1, format="iso", scale="utc")
    tmp_t2 = Time(time_2, format="iso", scale="utc")
    date_obs = time_1.replace(" ", "T")
    time_star = float(round(tmp_t1.mjd, 5))
    time_stop = float(round(tmp_t2.mjd, 5))
    time_stop = time_stop + (extCont[-1] / (24 * 60 * 60))
    time_stop = float("%0.5f" % time_stop) + 0.5
    stop_start = (time_stop - time_star) * (24 * 60 * 60)
    if ext_val > stop_start:
        print("")
        print(RED + "Exposure time cannot be greater than (MJD-END - MJD-OBS)" + ENDCOLOR)
        sys.exit()

    # --- ESO Programme identification ---
    print("\n" + CYAN + "Entry section" + ENDCOLOR)
    print(CYAN + ">>> " + ENDCOLOR + "Is this your ESO Programme Identification Code " + CYAN + str(prog_code) + ENDCOLOR + "?")
    if CONFIRM_PROG_CODE:
        print(GREEN + "ESO Programme Identification Code confirmed " + ENDCOLOR)
    else:
        prog_code = NEW_PROG_CODE
        print(GREEN + "ESO Programme Identification Code set to: " + prog_code + ENDCOLOR)

    # --- Bibliographic reference ---
    print("\n" + CYAN + ">>> " + ENDCOLOR + "Do you have Bibliographic reference for this data product?")
    print(CYAN + "NOTE" + ENDCOLOR + ": REFERENC must be the 19-digit bibliograp otherwise leave blank.")
    referen = BIBLIO_REF
    if referen:
        if len(referen) == 19 and referen[-1].isalpha():
            print(GREEN + "Bibcode keep it" + ENDCOLOR)
        else:
            print(RED + "Invalid Bibcode, entry left blank" + ENDCOLOR)
            referen = ""
    else:
        print(GREEN + "Bibcode entry left blank" + ENDCOLOR)

    # --- Compute additional header keywords ---
    crtoSofw = prihdr["ORIGIN"]
    refrq4ap = refreq / 1.0e9  # GHz
    WAVELMAX = float("%0.6f" % ((min(freqs) * u.GHz).to(u.nm, equivalencies=u.spectral()).value))
    WAVELMIN = float("%0.6f" % ((max(freqs) * u.GHz).to(u.nm, equivalencies=u.spectral()).value))
    lambda_c = (WAVELMIN + WAVELMAX) / 2.0
    bandwidth = WAVELMAX - WAVELMIN
    channelwidth = bandwidth / (specNr - 1)
    specres = lambda_c / channelwidth
    print("Spectral resolving power: ", specres)
    D = 7.8 * (800 / refrq4ap)
    D = float("%0.8f" % ((D * u.arcsec).to(u.degree).value))
    beam_size = 0.5 * (prihdr["BMAJ"] + prihdr["BMIN"]) * 3600.0
    print(CYAN + "Effective beam size from header is: " + str(beam_size) + ENDCOLOR)
    if NEW_BEAM_SIZE is not None:
        beam_size = NEW_BEAM_SIZE
        print(GREEN + "New beam size set to: " + str(beam_size) + ENDCOLOR)

    print(prihdr)

    # --- Delete unused keywords ---
    for key in ["BUNIT", "DATAMIN", "DATAMAX"]: #, "LINE"
        if key in prihdr:
            del prihdr[key]

    equinox = 2000
    c_frame = "FK5"
    if equinox == 2000 and c_frame == "FK5":
        coordRef = c_frame
    else:
        print(RED + "Coordinate frame error." + ENDCOLOR)
        sys.exit()

    # --- Update Primary HDU header ---
    print("-------------------------Primary HDU-------------------------")
    prihdr.set("PRODCATG", "SCIENCE.CUBE", "Data product category", 1)
    prihdr.set("ORIGIN", "APEX", "Facility")
    prihdr.set("TELESCOP", "APEX-12m", "Telescope name", after="ORIGIN")
    prihdr.set("INSTRUME", instrument, "Instrument name", after="TELESCOP")
    prihdr.set("FEBE1", febe, "APEX frontend/backend combination", after="INSTRUME")
    prihdr.set("OBJECT", source, "Target designation", after="FEBE1")
    prihdr.set("EXTEND", True, "Extensions are present")
    prihdr.set("RA", ra, "[deg] Derived from native frame configuration", after="OBJECT")
    prihdr.set("DEC", dec, "[deg] Derived from native frame configuration", after="RA")
    prihdr.set("EQUINOX", int(equinox), "Standard FK5 (years)", after="DEC")
    prihdr.set("RADESYS", coordRef, "Coordinate reference frame", after="EQUINOX")
    prihdr.set("TIMESYS", "TAI", "Time system for MJD", after="RADESYS")
    prihdr.set("EXPTIME", ext_val, "[s] Total integration time per exposure", after="TIMESYS")
    prihdr.set("TEXPTIME", ext_val, "[s] Total integration time of all exposures", after="EXPTIME")
    prihdr.set("MJD-OBS", time_star, "Start of observations (days)", after="TEXPTIME")
    prihdr.set("MJD-END", time_stop, "End of observations (days)", after="MJD-OBS")
    prihdr.set("DATE-OBS", date_obs, "start of observation ISO 8601 format", after="MJD-END")
    if prog_code == "MULTI":
        prihdr.set("PROG_ID", prog_code, after="TEXPTIME")
        for i, p in enumerate(proList, start=1):
            prihdr.set("PROGID" + str(i), p, "ESO programme identification", before="OBID1")
    else:
        prihdr.set("PROG_ID", prog_code, "ESO programme identification", after="TEXPTIME")
    prihdr.set("SPEC_RES", specres, "Average spectral resolving power", after="PROG_ID")
    prihdr.set("WAVELMIN", WAVELMIN, "[nm] Minimum wavelength", after="SPEC_RES")
    prihdr.set("WAVELMAX", WAVELMAX, "[nm] Maximum wavelength", after="WAVELMIN")
    prihdr.set("PROCSOFT", crtoSofw, "Reduction software/system", before="DATE")
    prihdr.set("OBSTECH", obs_tech, "Technique of observation")
    prihdr.set("MAPMODE", "OTF", "APEX map mode")
    prihdr.set("FLUXCAL", "ABSOLUTE", "Characterises the flux calibration", after="OBSTECH")
    prihdr.set("BNOISE", rms_med_Jy, "Median rms (Jy)")
    for i, f in enumerate(fileList, start=1):
        prihdr.set("PROV" + str(i), f.replace(" ", ":") + ".fits", "Original science file")
    prihdr.set("NCOMBINE", NrFilesBol, "Number of scans")
    prihdr.set("JYFACTOR", factor, "[Jy/K] The Jansky to Kelvin conversion factor")
    prihdr.set("JYUNIT", "[Jy/K]", "The Jansky to Kelvin conversion factor unit")
    prihdr.set("JYUNCERT", fac_un, "[Jy/K] The Jansky to Kelvin conversion uncertainty")
    prihdr.set("SKY_RES", beam_size, "effective beam size (arcsec)")
    prihdr.set("REFERENC", referen, "Bibliographic reference")

    # --- Create Science extension header ---
    exthdr = hdu_fluxcube.header
    exthdr.set("HDUCLASS", "ESO", "class name (ESO format)")
    exthdr.set("HDUDOC", "SDP", "ESO Science Data Products standard")
    exthdr.set("HDUVERS", "SDP version 7", "version number")
    exthdr.set("HDUCLAS1", "IMAGE", "data classification")
    exthdr.set("HDUCLAS2", "DATA", "this extension contains the science data")
    exthdr.set("ERRDATA", "STAT_EXT", "pointer to the error extension")
    exthdr["RA"] = (prihdr["RA"], "[deg] Spectroscopic target position (J2000.0)")
    exthdr["DEC"] = (prihdr["DEC"], "[deg] Spectroscopic target position (J2000.0)")
    exthdr["OBJECT"] = (prihdr["OBJECT"], "Target designation")
    exthdr["NAXIS"] = prihdr["NAXIS"]
    exthdr["NAXIS1"] = prihdr["NAXIS1"]
    exthdr["NAXIS2"] = prihdr["NAXIS2"]
    exthdr["NAXIS3"] = prihdr["NAXIS3"]
    exthdr.set("BUNIT", b_unit, "", after="OBJECT")
    exthdr["CTYPE1"] = prihdr["CTYPE1"]
    exthdr["CRVAL1"] = prihdr["CRVAL1"]
    exthdr["CRPIX1"] = prihdr["CRPIX1"]
    if "CUNIT1" in prihdr and prihdr["CUNIT1"]:
        exthdr["CUNIT1"] = prihdr["CUNIT1"]
    else:
        exthdr["CUNIT1"] = "deg"
    exthdr["CTYPE2"] = prihdr["CTYPE2"]
    exthdr["CRVAL2"] = prihdr["CRVAL2"]
    exthdr["CRPIX2"] = prihdr["CRPIX2"]
    if "CUNIT2" in prihdr and prihdr["CUNIT2"]:
        exthdr["CUNIT2"] = prihdr["CUNIT2"]
    else:
        exthdr["CUNIT2"] = "deg"
    exthdr["CTYPE3"] = prihdr["CTYPE3"]
    exthdr["CRVAL3"] = prihdr["CRVAL3"]
    exthdr["CRPIX3"] = prihdr["CRPIX3"]
    if "CD3_3" in prihdr:
        exthdr["CD3_3"] = prihdr["CD3_3"]
        exthdr["CD1_3"] = prihdr.get("CD1_3", 0.0)
        exthdr["CD2_3"] = prihdr.get("CD2_3", 0.0)
        exthdr["CD3_1"] = prihdr.get("CD3_1", 0.0)
        exthdr["CD3_2"] = prihdr.get("CD3_2", 0.0)
    else:
        if "CDELT3" in prihdr:
            exthdr.set("CD3_3", prihdr["CDELT3"], "Transformation matrix element spectral resolution", after="CRPIX3")
            exthdr.set("CD1_3", 0.0, "Transformation matrix element", after="CD3_3")
            exthdr.set("CD2_3", 0.0, "Transformation matrix element", after="CD1_3")
            exthdr.set("CD3_1", 0.0, "Transformation matrix element", after="CD2_3")
            exthdr.set("CD3_2", 0.0, "Transformation matrix element", after="CD3_1")
        else:
            print(RED + "No spectral resolution info found." + ENDCOLOR)
            sys.exit()
    if "CUNIT3" in prihdr and prihdr["CUNIT3"]:
        exthdr["CUNIT3"] = prihdr["CUNIT3"]
    else:
        exthdr["CUNIT3"] = "m/s"

    if (("CD1_1" in prihdr) and ("CD1_2" in prihdr) and ("CD2_1" in prihdr) and ("CD2_2" in prihdr)):
        print("CDi_j transformation matrix is present")
    else:
        if "CROTA2" in prihdr and ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
            print("Computing transformation matrix elements from CDELT and CROTA2")
            cd12 = abs(-1.0 * prihdr["CDELT2"] * math.sin(prihdr["CROTA2"]))
            cd21 = abs(prihdr["CDELT1"] * math.sin(prihdr["CROTA2"]))
            exthdr.set("CD1_1", prihdr["CDELT1"] * math.cos(prihdr["CROTA2"]), "Transformation matrix element")
            exthdr.set("CD1_2", cd12, "Transformation matrix element")
            exthdr.set("CD2_1", cd21, "Transformation matrix element")
            exthdr.set("CD2_2", prihdr["CDELT2"] * math.cos(prihdr["CROTA2"]), "Transformation matrix element")
        else:
            if ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
                print("No rotation found: creating simple transformation matrix")
                exthdr.set("CD1_1", prihdr["CDELT1"], "Transformation matrix element")
                exthdr.set("CD1_2", 0.0, "Transformation matrix element")
                exthdr.set("CD2_1", 0.0, "Transformation matrix element")
                exthdr.set("CD2_2", prihdr["CDELT2"], "Transformation matrix element")
            else:
                print(RED + "No transformation matrix available." + ENDCOLOR)
                sys.exit()

    # # --- Compute the spatial transformation matrix ---
    # if all(k in prihdr for k in ("PC1_1", "PC1_2", "PC2_1", "PC2_2")):
    #     # Convert from PC to CD using the CDELT values
    #     exthdr.set("CD1_1", prihdr["CDELT1"] * prihdr["PC1_1"], "Computed from PC matrix")
    #     exthdr.set("CD1_2", prihdr["CDELT1"] * prihdr["PC1_2"], "Computed from PC matrix")
    #     exthdr.set("CD2_1", prihdr["CDELT2"] * prihdr["PC2_1"], "Computed from PC matrix")
    #     exthdr.set("CD2_2", prihdr["CDELT2"] * prihdr["PC2_2"], "Computed from PC matrix")
    #     print("PC matrix converted to CD matrix")
    # elif (("CD1_1" in prihdr) and ("CD1_2" in prihdr) and ("CD2_1" in prihdr) and ("CD2_2" in prihdr)):
    #     print("CDi_j transformation matrix is present")
    # elif "CROTA2" in prihdr and ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
    #     print("Computing transformation matrix elements from CDELT and CROTA2")
    #     cd12 = abs(-1.0 * prihdr["CDELT2"] * math.sin(prihdr["CROTA2"]))
    #     cd21 = abs(prihdr["CDELT1"] * math.sin(prihdr["CROTA2"]))
    #     exthdr.set("CD1_1", prihdr["CDELT1"] * math.cos(prihdr["CROTA2"]), "Transformation matrix element")
    #     exthdr.set("CD1_2", cd12, "Transformation matrix element")
    #     exthdr.set("CD2_1", cd21, "Transformation matrix element")
    #     exthdr.set("CD2_2", prihdr["CDELT2"] * math.cos(prihdr["CROTA2"]), "Transformation matrix element")
    # elif ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
    #     print("No rotation found: creating simple transformation matrix")
    #     exthdr.set("CD1_1", prihdr["CDELT1"], "Transformation matrix element")
    #     exthdr.set("CD1_2", 0.0, "Transformation matrix element")
    #     exthdr.set("CD2_1", 0.0, "Transformation matrix element")
    #     exthdr.set("CD2_2", prihdr["CDELT2"], "Transformation matrix element")
    # else:
    #     print(RED + "No transformation matrix available." + ENDCOLOR)
    #     sys.exit()

    exthdr["SPECSYS"] = prihdr["SPECSYS"]
    exthdr["RESTFREQ"] = prihdr["RESTFREQ"]
    exthdr["VELO-LSR"] = prihdr["VELO-LSR"]

    try: 
        exthdr["IMAGFREQ"] = prihdr["IMAGFREQ"]
    except:
        print('IMAGFREQ not found in header')
        prihdr["IMAGFREQ"] = prihdr["RESTFREQ"]
        exthdr["IMAGFREQ"] = prihdr["IMAGFREQ"]

    # --- Create RMS extension header ---
    extrmshdr = hdu_rmscube.header
    extrmshdr.set("HDUCLASS", "ESO", "class name (ESO format)")
    extrmshdr.set("HDUDOC", "SDP", "ESO Science Data Products standard")
    extrmshdr.set("HDUVERS", "SDP version 7", "version number")
    extrmshdr.set("HDUCLAS1", "IMAGE", "data classification")
    extrmshdr.set("HDUCLAS2", "ERROR", "this extension contains the error data")
    extrmshdr.set("HDUCLAS3", "RMSE", "error type")
    extrmshdr.set("SCIDATA", "DATA_EXT", "pointer to the science extension")
    extrmshdr["RA"] = (prihdr["RA"], "[deg] Spectroscopic target position (J2000.0)")
    extrmshdr["DEC"] = (prihdr["DEC"], "[deg] Spectroscopic target position (J2000.0)")
    extrmshdr["OBJECT"] = (prihdr["OBJECT"], "Target designation")
    extrmshdr["NAXIS"] = prihdr["NAXIS"]
    extrmshdr["NAXIS1"] = prihdr["NAXIS1"]
    extrmshdr["NAXIS2"] = prihdr["NAXIS2"]
    extrmshdr["NAXIS3"] = prihdr["NAXIS3"]
    extrmshdr.set("BUNIT", b_unit_rms, "", after="OBJECT")
    extrmshdr["CTYPE1"] = prihdr["CTYPE1"]
    extrmshdr["CRVAL1"] = prihdr["CRVAL1"]
    extrmshdr["CRPIX1"] = prihdr["CRPIX1"]
    extrmshdr["CUNIT1"] = exthdr["CUNIT1"]
    extrmshdr["CTYPE2"] = prihdr["CTYPE2"]
    extrmshdr["CRVAL2"] = prihdr["CRVAL2"]
    extrmshdr["CRPIX2"] = prihdr["CRPIX2"]
    extrmshdr["CUNIT2"] = exthdr["CUNIT2"]
    extrmshdr["CTYPE3"] = prihdr["CTYPE3"]
    extrmshdr["CRVAL3"] = prihdr["CRVAL3"]
    extrmshdr["CRPIX3"] = prihdr["CRPIX3"]
    extrmshdr["CD3_3"] = exthdr["CD3_3"]
    extrmshdr["CD1_3"] = exthdr["CD1_3"]
    extrmshdr["CD2_3"] = exthdr["CD2_3"]
    extrmshdr["CD3_1"] = exthdr["CD3_1"]
    extrmshdr["CD3_2"] = exthdr["CD3_2"]
    extrmshdr["CUNIT3"] = exthdr["CUNIT3"]
    extrmshdr["CD1_1"] = exthdr["CD1_1"]
    extrmshdr["CD1_2"] = exthdr["CD1_2"]
    extrmshdr["CD2_1"] = exthdr["CD2_1"]
    extrmshdr["CD2_2"] = exthdr["CD2_2"]
    extrmshdr["SPECSYS"] = prihdr["SPECSYS"]
    extrmshdr["RESTFREQ"] = prihdr["RESTFREQ"]
    extrmshdr["VELO-LSR"] = prihdr["VELO-LSR"]
    extrmshdr["IMAGFREQ"] = prihdr["IMAGFREQ"]

    # --- Create white-light image header ---
    whitehdr = hdu_white.header
    whitehdr.set("PRODCATG", "ANCILLARY.IMAGE.WHITELIGHT")
    whitehdr.set("NAXIS", 2, "")
    whitehdr["RA"] = (prihdr["RA"], "[deg] Spectroscopic target position (J2000.0)")
    whitehdr["DEC"] = (prihdr["DEC"], "[deg] Spectroscopic target position (J2000.0)")
    whitehdr["OBJECT"] = (prihdr["OBJECT"], "Target designation")
    whitehdr.set("BUNIT", b_unit, "", after="OBJECT")
    whitehdr["CTYPE1"] = prihdr["CTYPE1"]
    whitehdr["CRVAL1"] = prihdr["CRVAL1"]
    whitehdr["CRPIX1"] = prihdr["CRPIX1"]
    whitehdr["CUNIT1"] = exthdr["CUNIT1"]
    whitehdr["CTYPE2"] = prihdr["CTYPE2"]
    whitehdr["CRVAL2"] = prihdr["CRVAL2"]
    whitehdr["CRPIX2"] = prihdr["CRPIX2"]
    whitehdr["CUNIT2"] = exthdr["CUNIT2"]
    whitehdr["CD1_1"] = exthdr["CD1_1"]
    whitehdr["CD1_2"] = exthdr["CD1_2"]
    whitehdr["CD2_1"] = exthdr["CD2_1"]
    whitehdr["CD2_2"] = exthdr["CD2_2"]

    orighdr.set("PRODCATG", "ANCILLARY.CUBE")

    # --- Delete some primary header keywords no longer needed ---
    for pattern in ["CTYPE*...", "CRPIX*...", "CRVAL*...", "CDELT*...", "CROTA*...", "CUNIT1*...", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "RESTFREQ", "VELO-LSR", "IMAGFREQ"]:
        if pattern in prihdr:
            del prihdr[pattern]

    print(prihdr)
    print("")
    print(exthdr)
    print("")
    print(extrmshdr)
    print("")
    print(whitehdr)

    # --- Virtual output file names and header associations ---
    # outFluname = source + "_3DcubePh3.fits"
    # outWhitename = source + "_3DcubePh3_whitelight.fits"
    outFluname = FILENAME_OUT
    outWhitename = FILENAME_WHITE_OUT

    prihdr.set("ASSON1", outWhitename)
    if attlmv == "y":
        outLMVname = source + "_3Dcube_GILDASformat.lmv"
        prihdr.set("ASSON2", outLMVname)
        prihdr.set("ASSOC2", "ANCILLARY.CUBE.LMV")
        md5 = hashlib.md5(open(lmvIn, "rb").read()).hexdigest()
        prihdr.set("ASSOM2", md5)

    print(GREEN + "New FITS files with your calibrated data product are created! \n"
          "They contain the data (flux and rms cube) in a HDU structure and with a header compliant with the Phase 3 requirements.\n"
          "In addition a white-light image and a FITS file containing the flux cube with the original HDU structure is created.\n"
          "The data should then be ready to be included in a Phase 3 submission." + ENDCOLOR)
    if attlmv == "y":
        print(GREEN + "An IRAM GILDAS format .lmv cube has been associated and will also have to be provided in the Phase 3 submission." + ENDCOLOR)
    print(GREEN + "Fingers crossed..." + ENDCOLOR)

    # --- Write output files ---
    hdulist = fits.HDUList([fits.PrimaryHDU(header=prihdr), hdu_fluxcube, hdu_rmscube])
    hdulist[0].header = prihdr
    hdulist[1].header = exthdr
    hdulist[2].header = extrmshdr
    hdu_white.header = whitehdr

    hdulist.writeto(outFluname, checksum=True, overwrite=True)
    print(outFluname + " created...")
    hdu_white.writeto(outWhitename, checksum=True, overwrite=True)
    print(outWhitename + " created...")
    if attlmv == "y":
        os.system("cp -i " + lmvIn + " " + outLMVname)
        print(outLMVname + " created (copy of original GILDAS .lmv format cube)")

if __name__ == "__main__":
    main()