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

# Colors for printing messages
RED = "\033[31m"       # For errors
GREEN = "\033[32;3m"   # For outputs
BLUE = "\033[34;1m"    # Request inputs
MAGENTA = "\033[35m"   # Header
CYAN = "\033[36m"      # Main words
ENDCOLOR = "\033[0m"

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
    "SUPERCAM": [2014, 53, 8],  # NEW: same value as for FLASH345 and HET345 also for SUPERCAM
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

    print("\nQuerying the ESO TAP service at %s" % ESO_TAP_OBS)
    query = (
        "SELECT dp_id, exposure, prog_id, object, dp_tech, instrument, ra, dec \n"
        "FROM dbo.raw \n"
        "WHERE dp_id LIKE 'APEXHET.%' \n"
        "AND (prog_id LIKE '092.F-9315%' OR prog_id LIKE '193.C-0584%') \n"
        "AND object LIKE 'G301%' \n"
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
    Get instrument-backend input from user.
    """
    try:
        print("- Enter Instrument-Backend used (e.g. HET230-XFFTS2, see http://archive.eso.org/wdb/wdb/eso/apex/query)")
        febe_usr = input(" >>> ")
    except Exception:
        print(GREEN + "FEBE needed, try once again" + ENDCOLOR)
        sys.exit()
    return febe_usr


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
    # velo = prihdr["VELO-LSR"]  # Not used further
    crfreq = refreq * (1 - crvel / v_light)
    freq0 = crfreq - (crpix1 - 1) * cdelf
    freqs = []
    print("Frecuency axis in construction")
    for i in range(specNr):
        f1 = (freq0 + i * cdelf) / 1.0e9  # Convert to GHz
        freqs.append(f1)
    return freqs, refreq, cdelf


##########################################
# Main script execution
##########################################
def main():
    # Header output
    print("##########################################")
    print(MAGENTA + "#APEX-Phase 3 Tool" + ENDCOLOR)
    print(MAGENTA + "#Creation of FITS spectral cubes" + ENDCOLOR)
    print(MAGENTA + "#compliant with ESO Phase 3 requirements" + ENDCOLOR)
    print(MAGENTA + "#ESO-APEX 12m telescope" + ENDCOLOR)
    print(MAGENTA + "#contact: pvenegas@eso.org, tstanke@mpe.mpg.de" + ENDCOLOR)
    print("##########################################\n")

    # Important information
    print(RED + "Important information:" + ENDCOLOR, "\n",
          "- All the RAW science data files (one per scan number) used to produce this data product",
          "should be present in the current folder (in ARCFILE format with extension .fits). Do not",
          "include calibration data.", "\n",
          "- RAW science files can be downloaded from http://archive.eso.org/wdb/wdb/eso/apex/form.",
          'Run the query to find your scan(s) and then make sure you select the "Mark RAW" column',
          "before requesting your data.\n")

    ##########################################
    # List current FITS files and open science file
    ##########################################
    print("\n" + CYAN + "Your current FITS files" + ENDCOLOR)
    os.system("ls *.fits")
    print("")

    listOfFiles = os.listdir(".")
    print("\n" + CYAN + "FILES" + ENDCOLOR)
    specIn = input("- Which one is your FITS file containing your data product (mandatory) : >> ")
    fits.info(specIn)
    specfile = fits.open(specIn)

    ##########################################
    # Open FLUX FITS file and extract key header/data info
    ##########################################
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
        ra = prihdr["RA"]
        dec = prihdr["DEC"]

        print("\nObject name : ", source, "; Flux unit: ", b_unit, "\n")
    except Exception:
        print(RED + "Wrong FLUX FITS file, try it once again." + ENDCOLOR)
        sys.exit()

    sci_wcs = WCS(specIn)
    specfile.close()

    ##########################################
    # Open and check RMS cube
    ##########################################
    rmsIn = specIn.replace(".fits", "_rmscube.fits")
    print("")
    if rmsIn in listOfFiles:
        answ = input(BLUE + f"{rmsIn} exists, use this as rms cube? (y/n) >>> " + ENDCOLOR).lower()
        if answ != "y":
            print(CYAN + "Your current FITS files" + ENDCOLOR)
            os.system("ls *.fits")
            print("")
            rmsIn = input(BLUE + "  --- Which one is the fits file containing the error cube (mandatory)? >> " + ENDCOLOR)
    else:
        print(CYAN + "Your current FITS files" + ENDCOLOR)
        os.system("ls *.fits")
        print("")
        rmsIn = input(BLUE + "  --- Which one is the fits file containing the error cube (mandatory)? >> " + ENDCOLOR)

    if rmsIn not in listOfFiles:
        print(RED + f"file '{rmsIn}' containing error cube does not exist \n Please provide a valid file containing the error cube" + ENDCOLOR)
        sys.exit()

    fits.info(rmsIn)
    rmsfile = fits.open(rmsIn)

    print(CYAN + "... checking rms cube ..." + ENDCOLOR)
    try:
        rmsdat = rmsfile[0].data
        rmshdr = rmsfile[0].header
        rmsNr = rmsdat.shape[0]
        b_unit_rms = str(rmshdr["BUNIT"])
        source_rms = rmshdr["OBJECT"]

        if ("RA" not in prihdr) and ("RA" in prihdr["CTYPE1"]):
            prihdr["RA"] = prihdr["CRVAL1"]
        if ("DEC" not in prihdr) and ("DEC" in prihdr["CTYPE2"]):
            prihdr["DEC"] = prihdr["CRVAL2"]
        ra = prihdr["RA"]
        dec = prihdr["DEC"]

        print("\nObject name : ", source_rms, "; Flux unit: ", b_unit_rms, "\n")
    except Exception:
        print(RED + "Wrong rms FITS file, try it once again." + ENDCOLOR)
        sys.exit()

    if rmshdr["NAXIS"] != 3:
        print(RED + "rms file seems not to be a cube... \n    ---   Please give a different file!" + ENDCOLOR)
        sys.exit()
    else:
        print(GREEN + "      ... rms file is a cube ..." + ENDCOLOR)

    if (rmshdr["NAXIS1"] != prihdr["NAXIS1"]) or (rmshdr["NAXIS2"] != prihdr["NAXIS2"]) or (rmshdr["NAXIS3"] != prihdr["NAXIS3"]):
        print(RED + "Dimension of science cube: " + ENDCOLOR, prihdr["NAXIS1"], prihdr["NAXIS2"], prihdr["NAXIS3"])
        print(RED + "Dimension of rms cube    : " + ENDCOLOR, rmshdr["NAXIS1"], rmshdr["NAXIS2"], rmshdr["NAXIS3"])
        print(RED + "rms and science cube dimensions disagree... \n    ---   Please check or give a different file!" + ENDCOLOR)
        sys.exit()
    else:
        print(GREEN + "      ... rms and science cube dimensions agree ..." + ENDCOLOR)

    if (rmshdr["CRVAL1"] != prihdr["CRVAL1"] or rmshdr["CRPIX1"] != prihdr["CRPIX1"] or rmshdr["CDELT1"] != prihdr["CDELT1"] or
        rmshdr["CRVAL2"] != prihdr["CRVAL2"] or rmshdr["CRPIX2"] != prihdr["CRPIX2"] or rmshdr["CDELT2"] != prihdr["CDELT2"]):
        print(RED + "X/Y coordinate definition of science cube: " + ENDCOLOR,
              prihdr["CRVAL1"], prihdr["CRPIX1"], prihdr["CDELT1"],
              "; ", prihdr["CRVAL2"], prihdr["CRPIX2"], prihdr["CDELT2"])
        print(RED + "X/Y coordinate definition of rms cube    : " + ENDCOLOR,
              rmshdr["CRVAL1"], rmshdr["CRPIX1"], rmshdr["CDELT1"],
              "; ", rmshdr["CRVAL2"], rmshdr["CRPIX2"], rmshdr["CDELT2"])
        print(RED + "rms and science cube X/Y coordinate definitions disagree..." + ENDCOLOR)
        answ = input(BLUE + "Please indicate whether the difference is insignificant, and the science cube coordinate system should be copied to the rms cube (y/n) >>> " + ENDCOLOR).lower()
        if answ == "y":
            print(CYAN + "      ... copying X/Y coordinate system from science cube to rms cube ..." + ENDCOLOR)
            for key in ("CRVAL1", "CRPIX1", "CDELT1", "CRVAL2", "CRPIX2", "CDELT2"):
                rmshdr[key] = prihdr[key]
        else:
            print(RED + "       ... please check input rms file or provide a different one." + ENDCOLOR)
            sys.exit()
    else:
        print(GREEN + "      ... rms and science cube X/Y coordinate definitions agree ..." + ENDCOLOR)

    if (rmshdr["CRVAL3"] != prihdr["CRVAL3"] or rmshdr["CRPIX3"] != prihdr["CRPIX3"] or rmshdr["CDELT3"] != prihdr["CDELT3"]):
        print(RED + "Spectral axis definition of science cube: " + ENDCOLOR,
              prihdr["CRVAL3"], prihdr["CRPIX3"], prihdr["CDELT3"])
        print(RED + "Spectral axis definition of rms cube    : " + ENDCOLOR,
              rmshdr["CRVAL3"], rmshdr["CRPIX3"], rmshdr["CDELT3"])
        print(RED + "rms and science cube spectral axis definitions disagree... \n    ---   Please check or give a different file!" + ENDCOLOR)
        sys.exit()
    else:
        print(GREEN + "      ... rms and science cube spectral axis definitions agree ..." + ENDCOLOR)

    rmsfile.close()

    ##########################################
    # Associate .lmv cube if available (Optional)
    ##########################################
    attlmv = ""
    print("\n\n" + CYAN +
          "In case the spectral line cube was produced using the IRAM GILDAS software," + ENDCOLOR)
    print(CYAN + "we here provide the possibility to specify a .lmv file to be associated with" + ENDCOLOR)
    print(CYAN + "the Phase 3 compliant fits cube (OPTIONAL)." + ENDCOLOR)
    print("")
    lmvIn = specIn.replace(".fits", ".lmv")
    if lmvIn in listOfFiles:
        answ = input(BLUE + f"{lmvIn} exists, associate this as IRAM/GILDAS lmv cube? (y/n) >>> " + ENDCOLOR).lower()
        if answ == "y":
            attlmv = "y"
        else:
            answ2 = input(BLUE + "Would you like to associate any other .lmv file as IRAM/GILDAS lmv cube? (y/n) >>> " + ENDCOLOR).lower()
            if answ2 != "y":
                print(GREEN + "No .lmv file will be associated, continuing" + ENDCOLOR)
                attlmv = "n"
            else:
                print(CYAN + "Your current .lmv files (in this directory)" + ENDCOLOR)
                os.system("ls *.lmv")
                print("")
                lmvIn = input(BLUE + "  --- Which is the file containing the .lmv cube? >> " + ENDCOLOR)
                attlmv = "y"
    else:
        print(CYAN + "Your current .lmv files (in this directory):" + ENDCOLOR)
        os.system("ls *.lmv")
        print("")
        answ = input(BLUE + "Would you like to associate any .lmv file as IRAM/GILDAS lmv cube? (y/n) >>> " + ENDCOLOR).lower()
        if answ != "y":
            print(GREEN + "No .lmv file will be associated, continuing" + ENDCOLOR)
            attlmv = "n"
        else:
            print(CYAN + "Your current .lmv files (in this directory)" + ENDCOLOR)
            os.system("ls *.lmv")
            print("")
            lmvIn = input(BLUE + "  --- Which is the file containing the .lmv cube? >> " + ENDCOLOR)
            attlmv = "y"

    if attlmv == "y":
        if (os.path.exists(lmvIn)) and (lmvIn.endswith("lmv")):
            print(CYAN + f"file {lmvIn} exists; a copy will be associated to the Phase 3 cube \n (but needs to be uploaded to the archive separately!)" + ENDCOLOR)
        else:
            print(RED + f"file '{lmvIn}' does not exist or has wrong extension; \n Please provide a valid file containing the .lmv cube" + ENDCOLOR)
            sys.exit()

    ##########################################
    # Read TAP file and extract associated info
    ##########################################
    fileList = []
    extCont = []
    proList = []
    sourceList_tmp = []
    dprtech_tmp = []
    instrum = []
    posra = []
    posdec = []

    queryfile = tapQuery_Het(source)
    print(CYAN + "If you need to modified the TAP file, the time is now!" + ENDCOLOR)
    print()
    input("Please press the Enter key to resume\n")
    print()

    # Read the TAP file
    fileq = f"{source}.tap"
    with open(fileq, "r") as queryfile:
        qfile = [line.split() for line in queryfile][1:-1]

    print(BLUE + "After pause your TAP file length is:" + ENDCOLOR)
    print("Length = " + str(len(qfile)) + " rows\n")

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
    sourceList = list(set(sourceList_tmp))
    dprtech = list(set(dprtech_tmp))
    obs_tech = dprtech[0]
    instrument = instrum[0]
    posraval = float(posra[0])
    posdecval = float(posdec[0])

    print(CYAN + "These files will be associated with the cube:\n " + "\n".join(fileList) + ENDCOLOR)

    NrFilesBol = len(fileList)
    if NrFilesBol >= 1:
        print("\n" + GREEN + "Totat <APEXHET.*.fits> files >> " + ENDCOLOR, NrFilesBol)
        print(*fileList, sep="\t")
        print("")
    else:
        print("\n" + RED + "ESO ARCFILE data <APEXHET.*.fits> does not exist at the current Path" + ENDCOLOR)
        print(RED + "http://archive.eso.org/" + ENDCOLOR)
        sys.exit()

    proList = list(set(proList))
    if len(proList) == 1:
        prog_code = proList[0]
    elif len(proList) > 1:
        prog_code = "MULTI"
    print(prog_code)

    ext_val = sum(extCont)
    if ext_val < 5:
        print(RED + "There is an integration time inconsistency. Your ARCFILEs related to this program may be corrupted, Check your data" + ENDCOLOR)
        sys.exit()

    listTime = []
    for tmp_t in fileList:
        fileTime = tmp_t[8:18] + " " + tmp_t[19:31]
        listTime.append(fileTime)
        print(tmp_t, fileTime)

    ##########################################
    # FEBE, instrument, and Jy/K conversion factor
    ##########################################
    # The FEBE value is hardcoded below
    febe = "SUPERCAM-SCBE"
    ind = febe.index("-")
    front = febe[:ind]
    fac_yr = dicc[front][0]
    factor = dicc[front][1]
    fac_un = dicc[front][2]
    print(obs_tech, febe, front, factor)

    print(CYAN + "APEX spectral line cubes may be ingested in the archive either on a \nbrightness temperature scale (K) or flux density (Jy) scale." + ENDCOLOR)
    print("")
    print(CYAN + "In any case, a Jy/K conversion factor has to be provided" + ENDCOLOR)
    print("")
    print(CYAN + "- To obtain the conversion factor (Jy/K), please refer to the factor that best fits the time range of your observations." + ENDCOLOR)
    print("http://www.apex-telescope.org/telescope/efficiency/index.php")
    print("")
    print(CYAN + f"Currently used Jy/K factor for {front}: {factor}+/-{fac_un}  (Year: {fac_yr})" + ENDCOLOR)
    print("")
    answ = input(BLUE + "- Do you want to provide a different Jy/K factor? (y/n) >>> " + ENDCOLOR).lower()
    if answ == "y":
        try:
            factor = float(input(BLUE + "- Please, provide the new Jy/K factor (e.g. 40): >>>  " + ENDCOLOR))
        except Exception:
            print(RED + "Jy/K conversion factor must be a number" + ENDCOLOR)
            sys.exit()
        try:
            fac_un = float(input(BLUE + "- Please, provide the uncertainty : >>>  " + ENDCOLOR))
        except Exception:
            print(RED + "Jy/K uncertainty must be a number" + ENDCOLOR)
            sys.exit()
        print(GREEN + "Value set to: " + ENDCOLOR, factor, "+/-", fac_un)

    ##########################################
    # Flux conversion handling
    ##########################################
    try:
        if b_unit[0] == "K":
            print(CYAN + "- Your data product is currently expressed in terms of an antenna temperature scale [K]." + ENDCOLOR)
            print("")
            if b_unit != "K":
                print(CYAN + f"The brightness temperature unit is given in {b_unit}, \npossibly no correction for beam efficiency has been applied?" + ENDCOLOR)
                print("")
                print(CYAN + f"Please indicate whether a beam efficiency correction is still needed (if not, only the unit will be changed from '{b_unit}' to 'K')." + ENDCOLOR)
                answ = input(BLUE + " --- Do you want to apply a beam efficiency correction? (y/n) >>> " + ENDCOLOR).lower()
                if answ == "y":
                    try:
                        beam_eff = float(input(BLUE + "Please enter the beam efficiency (must be between 0 and 1 (greater than 0)) >>>" + ENDCOLOR))
                    except Exception:
                        print(RED + "Beam efficiency must be a number (between 0 and 1))" + ENDCOLOR)
                        sys.exit()
                    if 0.0 < beam_eff <= 1.0:
                        print(CYAN + "        ... applying a beam efficiency of: " + str(beam_eff) + ENDCOLOR)
                        print(CYAN + "        ... and changing brightness unit from '" + b_unit + "' to 'K'." + ENDCOLOR)
                        scidat = scidat / beam_eff
                        rmsdat = rmsdat / beam_eff
                    else:
                        print(RED + "Beam efficiency must be between 0 and 1 (greater than 0 ...)" + ENDCOLOR)
                        sys.exit()
                    print(beam_eff)
                else:
                    print(CYAN + "        ... Changing brightness unit from '" + b_unit + "' to 'K' (without applying any further beam efficiency correction)" + ENDCOLOR)
                b_unit = "K"
                b_unit_rms = "K"
            else:
                print(CYAN + "The brightness temperature unit is given in " + b_unit + ", probably a correction for beam efficiency has been applied?" + ENDCOLOR)
                answ = input(BLUE + "Please indicate whether this is correct (y) or whether you would like to apply an additional correction (n) >>> " + ENDCOLOR).lower()
                if answ != "y":
                    try:
                        beam_eff = float(input(BLUE + "Please enter the beam efficiency (must be between 0 and 1 (greater than 0)) >>>" + ENDCOLOR))
                    except Exception:
                        print(RED + "Beam efficiency must be a number (between 0 and 1))" + ENDCOLOR)
                        sys.exit()
                    if 0.0 < beam_eff <= 1.0:
                        print(CYAN + "        ... applying a beam efficiency of: " + str(beam_eff) + ENDCOLOR)
                        scidat = scidat / beam_eff
                        rmsdat = rmsdat / beam_eff
                    else:
                        print(RED + "Beam efficiency must be between 0 and 1 (greater than 0 ...)" + ENDCOLOR)
                        sys.exit()
                else:
                    print(CYAN + "        ... no further beam efficiency correction applied ..." + ENDCOLOR)
            print("\n\n" + CYAN + "APEX spectral line cubes may be ingested in the archive either on a brightness temperature scale (K) or flux density (Jy) scale." + ENDCOLOR)
            answ = input(BLUE + "  --- Please indicate whether you want to convert the cube from 'K' to 'Jy' (y/n) >>> " + ENDCOLOR).lower()
            if answ == "y":
                flux = scidat * factor
                flux_rms = rmsdat * factor
                b_unit = "Jy"
                b_unit_rms = "Jy"
            else:
                flux = scidat
                flux_rms = rmsdat
        elif b_unit == "mJy":
            print(CYAN + "- Your data product looks like it is calibrated in terms of Flux density [mJy]" + ENDCOLOR)
            answ = input(BLUE + "  --- Do you whish to convert to a brightness temperature (K) scale? (y/n) >>> " + ENDCOLOR).lower()
            if answ == "y":
                flux = scidat / (1000.0 * factor)
                flux_rms = rmsdat / (1000.0 * factor)
                b_unit = "K"
                b_unit_rms = "K"
            else:
                flux = scidat
                flux_rms = rmsdat
        elif b_unit == "Jy":
            print(CYAN + "- Your data product looks like it is calibrated in terms of Flux density [Jy]" + ENDCOLOR)
            print("")
            answ = input(BLUE + "  --- Do you whish to convert to a brightness temperature (K) scale? (y/n) >>> " + ENDCOLOR).lower()
            if answ == "y":
                flux = scidat / factor
                flux_rms = rmsdat / factor
                b_unit = "K"
                b_unit_rms = "K"
            else:
                flux = scidat
                flux_rms = rmsdat
        else:
            print(RED + "Your data seems to be neither in a brightness temperature (K) nor flux density (mJy/Jy) scale" + ENDCOLOR)
            print(RED + "Please check the calibration and provide a cube in proper flux/brightness units." + ENDCOLOR)
            sys.exit()
    except Exception:
        if b_unit[0] == "K":
            print(RED + "000 Intensity seems to be on a brightness temperature scale. Please try it once again." + ENDCOLOR)
            sys.exit()
        elif (b_unit == "Jy") or (b_unit == "mJy"):
            print(RED + "000 Intensity seems to be on a flux density scale. Please try it once again." + ENDCOLOR)
            sys.exit()
        else:
            print(RED + "000 Please check intensity units (should be on a 'K', 'Jy', or 'mJy' scale; currently it seems to be on a '" + b_unit + "' scale) and/or try again..." + ENDCOLOR)
            sys.exit()

    def_val = "%0.2f" % np.nanmax(flux)
    rms_med = np.nanmedian(flux_rms)
    rms_med_Jy = rms_med * factor if b_unit_rms[0] == "K" else rms_med
    flux_white = np.nanmean(flux, axis=0)
    freqs, refreq, cdel = freqArray(prihdr, specNr)

    # Create HDUs for flux cube, rms cube, and white-light image
    hdu_fluxcube = fits.ImageHDU(data=flux, header=None, name="DATA_EXT")
    hdu_rmscube = fits.ImageHDU(data=flux_rms, header=None, name="STAT_EXT")
    hdu_white = fits.PrimaryHDU(flux_white)

    ##########################################
    # Observation timing validation
    ##########################################
    time_1 = listTime[0]
    time_2 = listTime[-1]
    tmp_t1 = Time(time_1, format="iso", scale="utc")
    tmp_t2 = Time(time_2, format="iso", scale="utc")
    date_obs = time_1.replace(" ", "T")
    time_star = float(round(tmp_t1.mjd, 5))
    time_stop = float(round(tmp_t2.mjd, 5))
    time_stop = float(time_stop) + (extCont[-1] / (24 * 60 * 60))
    time_stop = float("%0.5f" % time_stop) + 0.5
    stop_start = (time_stop - time_star) * (24 * 60 * 60)
    if ext_val > stop_start:
        print("")
        print(RED + "Exposure time cannot be greater than (MJD-END - MJD-OBS)" + ENDCOLOR)
        print(str(ext_val) + "[s] ExpTime > " + str(int(stop_start)) + "[s] observing Time!")
        sys.exit()

    ##########################################
    # ESO Programme validation
    ##########################################
    print("\n" + CYAN + "Entry section" + ENDCOLOR)
    print(CYAN + ">>> " + ENDCOLOR + "Is this your ESO Programme Identification Code " + CYAN + str(prog_code) + ENDCOLOR + "?")
    answ = input(" y / n >>>  ").lower()
    if answ == "y":
        print(GREEN + "ESO Programme Identification Code confirmed " + ENDCOLOR)
    elif answ == "n":
        proj_id = input("Enter the ESO Programme Identification Code / format:TP.C-NNNN(R) >>>  ")
        if ("." in proj_id) and ("-" in proj_id) and ("(" in proj_id) and (")" in proj_id):
            prog_code = proj_id
        else:
            print(RED + "The code does not meet the minimum format requirements - EXIT" + ENDCOLOR)
            sys.exit()
    else:
        print(GREEN + "ESO Programme Identification Code keep it" + ENDCOLOR)

    ##########################################
    # Bibliographic reference (Bibcode) validation
    ##########################################
    print("\n" + CYAN + ">>> " + ENDCOLOR + "Do you have Bibliographic reference for this data product?")
    print(CYAN + "NOTE" + ENDCOLOR + ": REFERENC must be the 19-digit bibliograp otherwise press ENTER.")
    referen = input("Bibliographic reference: >>>  ")
    ref_len = len(referen)
    try:
        if ref_len == 0:
            print(GREEN + "Bibcode entry left blank" + ENDCOLOR)
        elif ref_len > 0:
            if ref_len == 19 and referen[-1].isalpha():
                print(GREEN + "Bibcode keep it" + ENDCOLOR)
            else:
                print(RED + "That's invalid Bibcode, please try again." + ENDCOLOR)
                referen = input("Bibliographic reference (YYYYJJJJJVVVVMPPPPA): >>>  ")
                if len(referen) == 0:
                    print(GREEN + "Bibcode entry left blank" + ENDCOLOR)
                elif len(referen) == 19 and referen[-1].isalpha():
                    print(GREEN + "Bibcode keep it" + ENDCOLOR)
                else:
                    print(RED + "That's invalid Bibcode, Bibcode entry left blank" + ENDCOLOR)
                    referen = ""
    except Exception:
        pass

    ##########################################
    # Keywords and header modifications
    ##########################################
    crtDate = prihdr["DATE"]
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
    print(CYAN + "Effective beam size as obtained from the header of the input science cube is: " + str(beam_size) + ENDCOLOR)
    answ = input(BLUE + "   --- Would you like to enter a different beam size, e.g., to account for additional smoothing (y/n)? >>> " + ENDCOLOR).lower()
    if answ == "y":
        try:
            beam_size = float(input(BLUE + "Enter new beam size >>> " + ENDCOLOR))
        except Exception:
            print(RED + "Beam size must be a number, please try again!" + ENDCOLOR)
            sys.exit()
    print(prihdr)

    # Delete keywords that are not needed any more
    for key in ["BUNIT", "DATAMIN", "DATAMAX", "LINE"]:
        if key in prihdr:
            del prihdr[key]

    equinox = 2000
    c_frame = "FK5"
    if equinox == 2000 and c_frame == "FK5":
        coordRef = c_frame
    else:
        print(RED + "Check your coord. frame. At APEX is FK5, why your it does not?" + ENDCOLOR)
        sys.exit()

    # Primary HDU keywords
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

    # Science extension header
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
    if "CUNIT1" in prihdr:
        exthdr["CUNIT1"] = prihdr["CUNIT1"]
        if not prihdr["CUNIT1"]:
            if "CDELT1" in prihdr:
                pixelsize = prihdr["CDELT1"] * 3600.0
            elif "CD1_1" in prihdr:
                pixelsize = prihdr["CD1_1"] * 3600.0
            else:
                print(RED + "... neither CDELT1 nor CD1_1 are in the input fits header. Please check the WCS headers of the input file" + ENDCOLOR)
                sys.exit()
            print(CYAN + "Assuming that CDELT1 (or CD1_1) is in degrees, the pixel size is " + str(pixelsize) + " arcsec" + ENDCOLOR)
            if input(BLUE + "Unit of first WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> " + ENDCOLOR).lower() == "y":
                exthdr["CUNIT1"] = "deg"
            else:
                print(RED + "CUNITs for celestial coordinate system MUST be degrees." + ENDCOLOR)
                sys.exit()
        else:
            if prihdr["CUNIT1"] != "deg":
                print(RED + "CUNITs for celestial coordinate system MUST be degrees." + ENDCOLOR)
                sys.exit()
    else:
        if "CDELT1" in prihdr:
            pixelsize = prihdr["CDELT1"] * 3600.0
        elif "CD1_1" in prihdr:
            pixelsize = prihdr["CD1_1"] * 3600.0
        else:
            print(RED + "... neither CDELT1 nor CD1_1 are in the input fits header. Please check the WCS headers of the input file" + ENDCOLOR)
            sys.exit()
        print(CYAN + "Assuming that CDELT1 (or CD1_1) is in degrees, the pixel size is " + str(pixelsize) + " arcsec" + ENDCOLOR)
        if input(BLUE + "Unit of first WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> " + ENDCOLOR).lower() == "y":
            exthdr.set("CUNIT1", "deg", "", after="CRPIX1")
        else:
            print(RED + "CUNIT for celestial coordinate system MUST be degrees." + ENDCOLOR)
            sys.exit()

    exthdr["CTYPE2"] = prihdr["CTYPE2"]
    exthdr["CRVAL2"] = prihdr["CRVAL2"]
    exthdr["CRPIX2"] = prihdr["CRPIX2"]
    if "CUNIT2" in prihdr:
        exthdr["CUNIT2"] = prihdr["CUNIT2"]
        if not prihdr["CUNIT2"]:
            if "CDELT2" in prihdr:
                pixelsize = prihdr["CDELT2"] * 3600.0
            elif "CD2_2" in prihdr:
                pixelsize = prihdr["CD2_2"] * 3600.0
            else:
                print(RED + "... neither CDELT2 nor CD2_2 are in the input fits header. Please check the WCS headers of the input file" + ENDCOLOR)
                sys.exit()
            print(CYAN + "Assuming that CDELT2 (or CD2_2) is in degrees, the pixel size is " + str(pixelsize) + " arcsec" + ENDCOLOR)
            if input(BLUE + "Unit of second WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> " + ENDCOLOR).lower() == "y":
                exthdr["CUNIT2"] = "deg"
            else:
                print(RED + "CUNITs for celestial coordinate system MUST be degrees." + ENDCOLOR)
                sys.exit()
        else:
            if prihdr["CUNIT2"] != "deg":
                print(RED + "CUNITs for celestial coordinate system MUST be degrees." + ENDCOLOR)
                sys.exit()
    else:
        if "CDELT2" in prihdr:
            pixelsize = prihdr["CDELT2"] * 3600.0
        elif "CD2_2" in prihdr:
            pixelsize = prihdr["CD2_2"] * 3600.0
        else:
            print(RED + "... neither CDELT2 nor CD2_2 are in the input fits header. Please check the WCS headers of the input file" + ENDCOLOR)
            sys.exit()
        print(CYAN + "Assuming that CDELT2 (or CD2_2) is in degrees, the pixel size is " + str(pixelsize) + " arcsec" + ENDCOLOR)
        if input(BLUE + "Unit of second WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> " + ENDCOLOR).lower() == "y":
            exthdr.set("CUNIT2", "deg", "", after="CRPIX2")
        else:
            print(RED + "CUNITs for celestial coordinate system MUST be degrees." + ENDCOLOR)
            sys.exit()

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
            print(RED + "Can't find any information on spectral resolution (CDELT3 or CD3_3), this seems not to be a cube?" + ENDCOLOR)
            sys.exit()

    if "CUNIT3" in prihdr:
        exthdr["CUNIT3"] = prihdr["CUNIT3"]
    else:
        print("\n" + RED + "CUNIT3 is not given." + ENDCOLOR)
        print(CYAN + "Spectral axis header values are:")
        print(CYAN + "Spectral axis type: " + prihdr["CTYPE3"] + ENDCOLOR)
        print(CYAN + "Reference channel value: " + str(prihdr["CRVAL3"]) + ENDCOLOR)
        print(CYAN + "Channel width: " + str(prihdr["CDELT3"]) + ENDCOLOR)
        if prihdr["CTYPE3"] == "VRAD":
            print(CYAN + "Frequency axis is in radial velocity units." + ENDCOLOR)
            answ = input(BLUE + "Unit of frequency axis will be set to 'm/s'. Is that ok? (y/n) >>> " + ENDCOLOR).lower()
            funits = "m/s" if answ == "y" else input(BLUE + "Please enter units for the frequency axis (e.g., m/s, km/s) >>> " + ENDCOLOR)
        else:
            funits = input(BLUE + "Please enter units for the frequency axis (e.g., m/s, km/s, GHz) >>> " + ENDCOLOR)
        print("Frequency axis units will be set to: " + funits)
        exthdr.set("CUNIT3", funits, "", after="CRPIX3")

    if (("CD1_1" in prihdr) and ("CD1_2" in prihdr) and ("CD2_1" in prihdr) and ("CD2_2" in prihdr)):
        print("CDi_j transformation matrix is present")
    else:
        print("CDi_j transformation matrix is not present in input header.")
        if "CROTA2" in prihdr:
            if ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
                print("Computing CDi_j transformation matrix elements from CDELTi and CROTA2 keywords")
                cd12 = -1.0 * prihdr["CDELT2"] * math.sin(prihdr["CROTA2"])
                cd21 = prihdr["CDELT1"] * math.sin(prihdr["CROTA2"])
                cd12 = abs(cd12)
                cd21 = abs(cd21)
                exthdr.set("CD1_1", prihdr["CDELT1"] * math.cos(prihdr["CROTA2"]), "Transformation matrix element")
                exthdr.set("CD1_2", cd12, "Transformation matrix element")
                exthdr.set("CD2_1", cd21, "Transformation matrix element")
                exthdr.set("CD2_2", prihdr["CDELT2"] * math.cos(prihdr["CROTA2"]), "Transformation matrix element")
            else:
                print(RED + "Neither transformation matrix CDi_j nor CDELT1/CDELT2 are given" + ENDCOLOR)
                sys.exit()
        else:
            if ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
                print("Neither transformation matrix nor CROTA2 keyword found.")
                print("Creating transformation matrix assuming that there is no rotation to the coordinate system")
                exthdr.set("CD1_1", prihdr["CDELT1"], "Transformation matrix element")
                exthdr.set("CD1_2", 0.0, "Transformation matrix element")
                exthdr.set("CD2_1", 0.0, "Transformation matrix element")
                exthdr.set("CD2_2", prihdr["CDELT2"], "Transformation matrix element")
            else:
                print(RED + "Neither transformation matrix CDi_j nor CDELT1/CDELT2 are given" + ENDCOLOR)
                sys.exit()

    exthdr["SPECSYS"] = prihdr["SPECSYS"]
    exthdr["RESTFREQ"] = prihdr["RESTFREQ"]
    exthdr["VELO-LSR"] = prihdr["VELO-LSR"]
    exthdr["IMAGFREQ"] = prihdr["IMAGFREQ"]

    # RMS extension header
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

    # White-light image header
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

    orighdr.set("PRODCATG", "ANCILLARY.CUBE")

    # Delete primary header keywords that were moved
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

    ##########################################
    # Virtual output file names and header association
    ##########################################
    outFluname = source + "_3DcubePh3.fits"
    outWhitename = source + "_3DcubePh3_whitelight.fits"
    prihdr.set("ASSON1", outWhitename)
    if attlmv == "y":
        outLMVname = source + "_3Dcube_GILDASformat.lmv"
        prihdr.set("ASSON2", outLMVname)
        prihdr.set("ASSOC2", "ANCILLARY.CUBE.LMV")
        md5 = hashlib.md5(open(lmvIn, "rb").read()).hexdigest()
        prihdr.set("ASSOM2", md5)

    print(GREEN + "New FITS files with your calibrated data product are created!. \nThey contain the data (flux and rms cube) in a HDU structure and with a header compliant with the Phase 3 requirements.\nIn addition a white-light image and a fits file containing the flux cube with the original HDU structure is created.\nThe data should then be ready to be included in a Phase 3 submission." + ENDCOLOR)
    if attlmv == "y":
        print(GREEN + "An IRAM GILDAS format .lmv cube has been associated and will also have to be provided in the Phase 3 submission." + ENDCOLOR)
    print(GREEN + "Fingers crossed..." + ENDCOLOR)

    # Files creation
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