#!/usr/bin/env python
##########################################
# 2020-07-13		Creation Version 1.
# 2020-08-28		Pipeline update. Jy/K factor mandatory
# 2020-09-19		Keywords ASSOCi, ASSONi added in this pipeline.
# 2021-08-10             Removing keyword OBID1, which is not is not needed and makes the Ph3 fail verification
# 2021-08-13             Modifying keywords, adding capability to read tap files instead of raw fits files
##########################################
# BLACK	 = "\033[30m"
RED = "\033[31m"  # For errors
GREEN = "\033[32;3m"  # For outputs
# YELLOW   = "\033[33;1m"	# Not in use
BLUE = "\033[34;1m"  # Request inputs
MAGENTA = "\033[35m"  # Header
CYAN = "\033[36m"  # Main  words
# WHITE    = "\033[37m"
ENDCOLOR = "\033[0m"
##########################################
print("##########################################")
print(MAGENTA + "#APEX-Phase 3 Tool" + ENDCOLOR)
print(MAGENTA + "#Creation of FITS spectral cubes" + ENDCOLOR)
print(MAGENTA + "#compliant with ESO Phase 3 requirements" + ENDCOLOR)
print(MAGENTA + "#ESO-APEX 12m telescope" + ENDCOLOR)
print(MAGENTA + "#contact: pvenegas@eso.org, tstanke@mpe.mpg.de" + ENDCOLOR)
print("##########################################")
print("")
##########################################
# Libraries:
import numpy as np
import math
import os
import sys
import pathlib
import re
import datetime
import copy
import hashlib
from astropy.time import Time
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
from astropy.wcs import Wcsprm

from gadgets import functions

v_light = 299792458.0  # Speed of light (m/s)
##########################################

###Information
print(
    RED + "Important information:" + ENDCOLOR,
    "\n",
    "- All the RAW science data files (one per scan number) used to produce this data product",
    "should be present in the current folder (in ARCFILE format with extension .fits). Do not",
    "include calibration data.",
    "\n",
    "- RAW science files can be downloaded from http://archive.eso.org/wdb/wdb/eso/apex/form.",
    'Run the query to find your scan(s) and then make sure you select the "Mark RAW" column',
    "before requesting your data.",
)

# ----------------------------------------------------------------------------------------
################################
# Diccionario
################################
# Instrum - year - Jy/K - uncertainty
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
    "PI230": [2019, 35, 3],
}


dicc = {
    "HET230": [2016, 40, 6],
    "PI230": [2017, 45, 7],
    "HET345": [2016, 53, 8],
    # NEW: use same value as for FLASH345 and HET345 also for SUPERCAM
    "SUPERCAM": [2014, 53, 8],
    "FLASH345": [2018, 53, 4],
    "HET460": [2017, 72, 11],
    "GARD180": [2017, 40, 6],
    "NOVA660": [2017, 110, 18],
    "SEPIA180": [2019, 35, 3],
    "SEPIA660": [2019, 68, 6],
}

################################
# Definitions:
################################


########
def tapQuery_Het(source):
    # Version: ALU 2021-08-11

    import sys
    import math
    from pyvo.dal import tap

    from astropy.coordinates import SkyCoord
    from astropy.units import Quantity
    from tabulate import tabulate

    ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
    tapobs = tap.TAPService(ESO_TAP_OBS)

    target = source

    print()
    print("Querying the ESO TAP service at %s" % (ESO_TAP_OBS))

    # --------------------------------------------------
    # The actual position of the selected target
    # is queried by the from_name() function,
    # which queries the CDS SESAME service
    # (http://cdsweb.u-strasbg.fr/cgi-bin/Sesame).
    # --------------------------------------------------

    # query = """SELECT dp_id, exposure, prog_id, object, dp_tech, instrument, ra, dec 
	# from dbo.raw
	# where dp_id like 'APEXHET.%%'
  	# 	and object = '%s'
  	# 	and dp_cat='SCIENCE'""" % (target)
    
    query = """SELECT dp_id, exposure, prog_id, object, dp_tech, instrument, ra, dec 
    FROM dbo.raw
    WHERE dp_id LIKE 'APEXHET.%'
    AND (prog_id LIKE '092.F-9315%' OR prog_id LIKE '193.C-0584%')
    AND object LIKE 'G301%'
    AND dp_cat = 'SCIENCE'
    """ 
           
    print()
    print("Query:")
    print(query)

    res = tapobs.search(query=query, maxrec=1000)

    print()
    print(res.to_table())
    print()

    print()
    print(
        "A total of "
        + str(len(res.to_table()))
        + " records were found matching the provided criteria."
    )

    table = res.to_table()

    filename = "%s.tap" % target
    print()
    print("Results has been written to: %s" % filename)

    with open(filename, "w") as f:
        print(tabulate(table), file=f)

    f.close()
    return table


################################
def get_frontBack():
    ########
    # Description
    try:
        print(
            "- Enter Instrument-Backend used (e.g. HET230-XFFTS2, see http://archive.eso.org/wdb/wdb/eso/apex/query)"
        )
        febe_usr = str(input(" >>> "))

    except:
        print(GREEN + "FEBE needed, try once again" + ENDCOLOR)
        sys.exit()

    return febe_usr


################################
def freqArray(prihdr, specNr):
    ## Frequency array from the WCS
    cdelv = prihdr["CDELT3"]
    crpix1 = prihdr["CRPIX3"]
    crvel = prihdr["CRVAL3"]
    refreq = prihdr["RESTFREQ"]
    cdelf = -1.0 * refreq * cdelv / v_light
    velo = prihdr["VELO-LSR"]
    crfreq = refreq * (1 - crvel / v_light)

    freq0 = crfreq - (crpix1 - 1) * cdelf
    freqs = []
    print("Frecuency axis in construction")
    for i in range(specNr):
        f1 = (freq0 + i * cdelf) / 1.0e9  # == 10**9, GHz
        # print (f1)
        freqs.append(f1)
    return (freqs, refreq, cdelf)


# ----------------------------------------------------------------------------------------
################################
# Calling files in the current folder
################################
print("")
print(CYAN + "Your current FITS files" + ENDCOLOR)
os.system("ls *.fits")
print("")

listOfFiles = os.listdir(".")
# pattern = "APEX-*.fits"

################################
###OPEN Main files.
################################
print("")
print(CYAN + "FILES" + ENDCOLOR)
specIn = str(
    input(
        "- Which one is your FITS file containing your data product (mandatory) : >> "
    )
)
# specIn = 'cube-L1622-dsp-rsc-fx-fbl-bl-dsp-bl-dsp-bl-dsp-bl-dsp-fx-bl-bef-rsa-bl.fits'
fits.info(specIn)
specfile = fits.open(specIn)

################################
# OPEN FLUX Fits files
################################
try:
    scidat = specfile[0].data
    prihdr = specfile[0].header
    orighdr = copy.copy(prihdr)
    print(prihdr[0])
    specNr = scidat.shape[0]  # NEW this makes more sense for a cube
    b_unit = str(prihdr["BUNIT"])
    source = prihdr["OBJECT"]

    # extract ra and dec from the input fits file
    if ("RA" not in prihdr) & ("RA" in prihdr["CTYPE1"]):
        prihdr["RA"] = prihdr["CRVAL1"]
        
    if ("DEC" not in prihdr) & ("DEC" in prihdr["CTYPE2"]):
        prihdr["DEC"] = prihdr["CRVAL2"]
        
    ra = prihdr["RA"]
    dec = prihdr["DEC"]

    print("")
    print("Object name : ", source, "; Flux unit: ", b_unit)
    print("")
    cType1 = prihdr["CTYPE1"]

except:
    print(RED + "Wrong FLUX FITS file, try it once again." + ENDCOLOR)
    sys.exit()

sci_wcs = WCS(specIn)
# sci_wcs.wcs.print_contents()
specfile.close()

################################
###OPEN rms cube.
################################

rmsIn = copy.copy(specIn)
rmsIn = rmsIn.replace(".fits", "_rmscube.fits")
print("")
if rmsIn in listOfFiles:
    answ = str(
        input(BLUE + rmsIn + " exists, use this as rms cube? (y/n) >>> " + ENDCOLOR)
    )
    answ = answ.lower()
    if answ != "y":
        print(CYAN + "Your current FITS files" + ENDCOLOR)
        os.system("ls *.fits")
        print("")
        rmsIn = str(
            input(
                BLUE
                + "  --- Which one is the fits file containing the error cube (mandatory)? >> "
                + ENDCOLOR
            )
        )
else:
    print(CYAN + "Your current FITS files" + ENDCOLOR)
    os.system("ls *.fits")
    print("")
    rmsIn = str(
        input(
            BLUE
            + "  --- Which one is the fits file containing the error cube (mandatory)? >> "
            + ENDCOLOR
        )
    )

if rmsIn in listOfFiles:
    print(CYAN + "file " + rmsIn + " containing error cube exists" + ENDCOLOR)
else:
    print(
        RED
        + "file '"
        + rmsIn
        + "' containing error cube does not exist \n Please provide a valid file containing the error cube"
    )
    sys.exit()

fits.info(rmsIn)
rmsfile = fits.open(rmsIn)

### Check consistency of rms cube

print(CYAN + "... checking rms cube ..." + ENDCOLOR)

try:
    rmsdat = rmsfile[0].data
    rmshdr = rmsfile[0].header
    rmsNr = rmsdat.shape[0]
    b_unit_rms = str(rmshdr["BUNIT"])
    source_rms = rmshdr["OBJECT"]
    # extract ra and dec from the input fits file
    if ("RA" not in prihdr) & ("RA" in prihdr["CTYPE1"]):
        prihdr["RA"] = prihdr["CRVAL1"]
        
    if ("DEC" not in prihdr) & ("DEC" in prihdr["CTYPE2"]):
        prihdr["DEC"] = prihdr["CRVAL2"]
        
    ra = prihdr["RA"]
    dec = prihdr["DEC"]

    print("")
    print("Object name : ", source_rms, "; Flux unit: ", b_unit_rms)
    print("")
    cType1_rms = rmshdr["CTYPE1"]

except:
    print(RED + "Wrong rms FITS file, try it once again." + ENDCOLOR)
    sys.exit()


if rmshdr["NAXIS"] != 3:
    print(
        RED
        + "rms file seems not to be a cube... \n    ---   Please give a different file!"
        + ENDCOLOR
    )
    sys.exit()
else:
    print(GREEN + "      ... rms file is a cube ..." + ENDCOLOR)

if (
    (rmshdr["NAXIS1"] != prihdr["NAXIS1"])
    or (rmshdr["NAXIS2"] != prihdr["NAXIS2"])
    or (rmshdr["NAXIS3"] != prihdr["NAXIS3"])
):
    print(
        RED + "Dimension of science cube: " + ENDCOLOR,
        prihdr["NAXIS1"],
        prihdr["NAXIS2"],
        prihdr["NAXIS3"],
    )
    print(
        RED + "Dimension of rms cube    : " + ENDCOLOR,
        rmshdr["NAXIS1"],
        rmshdr["NAXIS2"],
        rmshdr["NAXIS3"],
    )
    print(
        RED
        + "rms and science cube dimensions disagree... \n    ---   Please check or give a different file!"
        + ENDCOLOR
    )
    sys.exit()
else:
    print(GREEN + "      ... rms and science cube dimensions agree ..." + ENDCOLOR)

if (
    (rmshdr["CRVAL1"] != prihdr["CRVAL1"])
    or (rmshdr["CRPIX1"] != prihdr["CRPIX1"])
    or (rmshdr["CDELT1"] != prihdr["CDELT1"])
    or (rmshdr["CRVAL2"] != prihdr["CRVAL2"])
    or (rmshdr["CRPIX2"] != prihdr["CRPIX2"])
    or (rmshdr["CDELT2"] != prihdr["CDELT2"])
):
    print(
        RED + "X/Y coordinate definition of science cube: " + ENDCOLOR,
        prihdr["CRVAL1"],
        prihdr["CRPIX1"],
        prihdr["CDELT1"],
        "; ",
        prihdr["CRVAL2"],
        prihdr["CRPIX2"],
        prihdr["CDELT2"],
    )
    print(
        RED + "X/Y coordinate definition of rms cube    : " + ENDCOLOR,
        rmshdr["CRVAL1"],
        rmshdr["CRPIX1"],
        rmshdr["CDELT1"],
        "; ",
        rmshdr["CRVAL2"],
        rmshdr["CRPIX2"],
        rmshdr["CDELT2"],
    )
    print(RED + "Central coordinate of science cube: " + ENDCOLOR, ra, dec)
    print(RED + "Central coordinate of rms cube    : " + ENDCOLOR, ra_rms, dec_rms)
    print(
        RED + "rms and science cube X/Y coordinate definitions disagree..." + ENDCOLOR
    )
    answ = str(
        input(
            BLUE
            + "Please indicate whether the difference is insignificant, and the science cube coordinate system should be copied to the rms cube (y/n) >>> "
            + ENDCOLOR
        )
    )
    answ = answ.lower()
    if answ == "y":
        print(
            CYAN
            + "      ... copying X/Y coordinate system from science cube to rms cube ..."
            + ENDCOLOR
        )
        rmshdr["CRVAL1"] = prihdr["CRVAL1"]
        rmshdr["CRPIX1"] = prihdr["CRPIX1"]
        rmshdr["CDELT1"] = prihdr["CDELT1"]
        rmshdr["CRVAL2"] = prihdr["CRVAL2"]
        rmshdr["CRPIX2"] = prihdr["CRPIX2"]
        rmshdr["CDELT2"] = prihdr["CDELT2"]
    else:
        print(
            RED
            + "       ... please check input rms file or provide a different one."
            + ENDCOLOR
        )
        sys.exit()
else:
    print(
        GREEN
        + "      ... rms and science cube X/Y coordinate definitions agree ..."
        + ENDCOLOR
    )

if (
    (rmshdr["CRVAL3"] != prihdr["CRVAL3"])
    or (rmshdr["CRPIX3"] != prihdr["CRPIX3"])
    or (rmshdr["CDELT3"] != prihdr["CDELT3"])
):
    print(
        RED + "Spectral axis definition of science cube: " + ENDCOLOR,
        prihdr["CRVAL3"],
        prihdr["CRPIX3"],
        prihdr["CDELT3"],
    )
    print(
        RED + "Spectral axis definition of rms cube    : " + ENDCOLOR,
        rmshdr["CRVAL3"],
        rmshdr["CRPIX3"],
        rmshdr["CDELT3"],
    )
    print(
        RED
        + "rms and science cube spectral axis definitions disagree... \n    ---   Please check or give a different file!"
        + ENDCOLOR
    )
    sys.exit()
else:
    print(
        GREEN
        + "      ... rms and science cube spectral axis definitions agree ..."
        + ENDCOLOR
    )


rmsfile.close()

# sys.exit()
# ----------------------------------------------------------------------------------------
################################
### Check whether there is a corresponding .lmv cube that should be associated
### (to provide a fully IRAM/GILDAS compatible dataset)
################################

attlmv = ""
print("")
print("")
print(
    CYAN
    + "In case the spectral line cube was produced using the IRAM GILDAS software,"
    + ENDCOLOR
)
print(
    CYAN
    + "we here provide the possibility to specify a .lmv file to be associated with"
    + ENDCOLOR
)
print(CYAN + "the Phase 3 compliant fits cube (OPTIONAL)." + ENDCOLOR)
print("")
lmvIn = copy.copy(specIn)
lmvIn = lmvIn.replace(".fits", ".lmv")
if lmvIn in listOfFiles:
    answ = str(
        input(
            BLUE
            + lmvIn
            + " exists, associate this as IRAM/GILDAS lmv cube? (y/n) >>> "
            + ENDCOLOR
        )
    )
    answ = answ.lower()
    if answ == "y":
        attlmv = "y"
    else:
        answ2 = str(
            input(
                BLUE
                + "Would you like to associate any other .lmv file as IRAM/GILDAS lmv cube? (y/n) >>> "
                + ENDCOLOR
            )
        )
        answ2 = answ2.lower()
        if answ2 != "y":
            print(GREEN + "No .lmv file will be associated, continuing" + ENDCOLOR)
            attlmv = "n"
        else:
            print(CYAN + "Your current .lmv files (in this directory)" + ENDCOLOR)
            os.system("ls *.lmv")
            print("")
            lmvIn = str(
                input(
                    BLUE
                    + "  --- Which is the file containing the .lmv cube? >> "
                    + ENDCOLOR
                )
            )
            attlmv = "y"
else:
    print(CYAN + "Your current .lmv files (in this directory):" + ENDCOLOR)
    os.system("ls *.lmv")
    print("")
    answ = str(
        input(
            BLUE
            + "Would you like to associate any other .lmv file as IRAM/GILDAS lmv cube? (y/n) >>> "
            + ENDCOLOR
        )
    )
    answ = answ.lower()
    if answ != "y":
        print(GREEN + "No .lmv file will be associated, continuing" + ENDCOLOR)
        attlmv = "n"
    else:
        print(CYAN + "Your current .lmv files (in this directory)" + ENDCOLOR)
        os.system("ls *.lmv")
        print("")
        lmvIn = str(
            input(
                BLUE
                + "  --- Which is the file containing the .lmv cube? >> "
                + ENDCOLOR
            )
        )
        attlmv = "y"

if attlmv == "y":
    if (os.path.exists(lmvIn)) and (lmvIn[-3:] == "lmv"):
        print(
            CYAN
            + "file "
            + lmvIn
            + " exists; a copy will be associated to the Phase 3 cube \n (but needs to be uploaded to the archive separately!)"
            + ENDCOLOR
        )
    else:
        print(
            RED
            + "file '"
            + lmvIn
            + "' does not exist or has wrong extension; \n Please provide a valid file containing the .lmv cube"
        )
        sys.exit()


################################
# OPEN APEXHET* files.
################################
fileList_tmp = []
fileList = []
sourceList_tmp = []
sourceList = []
extCont = []
listName = []
listTime = []
proList = []

### TAP NEW
dprtech_tmp = []
dprtech = []
instrum = []
posra = []
posdec = []

queryfile = tapQuery_Het(source)

print(CYAN + "If you need to modified the TAP file, the time is now!" + ENDCOLOR)
print()
input("Please press the Enter key to resume")
print()

####################
##Patch. To read the TAP file and not the table.
####################
# Open the TAP script.
fileq = "%s.tap" % source
queryfile = open(fileq, "r")

# Reading the TAP FILE
qfile = []
count = 0
for line in queryfile:
    # print(line, end="")
    qfile.append(line.split())
    count += 1

qfile = qfile[0:][1:-1]

count = 0
for line in qfile:
    # print (line)
    count += 1

print(BLUE + "After pause your TAP file length is:" + ENDCOLOR)
print("Length = " + str(count) + " rows")
print("")
queryfile.close()

for i in range(len(qfile)):
    fileList.append((qfile[i][0]))  # List of Fits files
    extCont.append(float(qfile[i][1]))  # Exposure Time
    proList.append(str(qfile[i][2]))  # Project Id.
    sourceList_tmp.append(str(qfile[i][3]))
    dprtech_tmp.append(str(qfile[i][4]))  # ESO DPR TECH / Scan technique
    instrum.append(str(qfile[i][5]))
    posra.append((qfile[i][6]))
    posdec.append((qfile[i][7]))

fileList.sort()
sourceList = list(set(sourceList_tmp))
dprtech = list(
    set(dprtech_tmp)
)  # TODO: add check that only one dprtech was used for the observations, anything else probably wouldn't make sense and point at an error in the tap query result
obs_tech = str(dprtech[0])
instrument = str(instrum[0])
posraval = float(posra[0])
posdecval = float(posdec[0])
### END NEW TAP
# sys.exit()

##### CONSTRUCTING LIST OF APEXHET FITS FILES PRESENT IN THE WORKING DIRECTORY,
##### CHECKING FOR CONSISTENCY, EXTRACTING EXPOSURE TIME, PROJECT IDS, ETC;
##### THIS IS NOT NEEDED ANYMORE, PROBABLY
## definition path
## currentDirectory = pathlib.Path(path)
# currentDirectory = pathlib.Path('.')

## define the pattern
# currentPattern = "APEXHET*.fits"

# for currentFile in currentDirectory.glob(currentPattern):
#    print(currentFile)
#    fileList_tmp.append(str(currentFile))
# fileList_tmp.sort()

# if len(fileList_tmp) < 1:
#        print(RED+"No APEXHET raw fits files found."+ENDCOLOR)
#        print(RED+"Please make sure that all APEXHET fits files pertaining to the cube are present in the current directory!"+ENDCOLOR)
#        print(RED+"(Download them from the archive and/or copy them to the current directory!)"+ENDCOLOR)
#        sys.exit()

### Figure out whether the source name in the cube equals the source name(s) in the APEXHET files
### or whether maps with different field names went into the final cube

# for apexhet in fileList_tmp:
# with fits.open(apexhet) as hdul:
#                hetinf = hdul[0].header
#                sub_sou = str(hetinf['OBJECT'])
#                sourceList_tmp.append(sub_sou)

# sourceList = list(set(sourceList_tmp))

#### THE FOLLOWING SECTION WAS USED TO CHECK WHETHER OR NOT THE APEXHET FITS FILES PRESENT IN THE
#### WORKING DIRECTORY DO OR DO NOT CORRESPOND TO THE FITS FILES THAT WERE USED FOR THE REDUCED CUBE
#### THIS NEEDS (?) TO BE ADAPTED TO WORK ON THE TAP QUERY INSTEAD OF THE LIST OF FITS FILES
#### (UNLESS WE CAN ASSUME THAT THE .TAP FILE HAS BEEN CHECKED FOR CONSISTENCY ALREADY)
# if len(sourceList) == 1:
#        if source.upper() == sourceList[0].upper():
#                print(GREEN+"All files share the same source name, continuing"+ENDCOLOR)
#                for apexhet in fileList_tmp:
# with fits.open(apexhet) as hdul:
#                                hetinf = hdul[0].header
#                                fileList.append(apexhet)
#                                extCont.append(int(hetinf['EXPTIME']))
#                                proList.append(str(hetinf['HIERARCH ESO OBS PROG ID']))
#        else:
#                print(RED+"Source name of cube does not correspond to source name given in APEXHET fits files."+ENDCOLOR)
#                print(RED+"Source name as given in cube: "+source+ENDCOLOR)
#                print(RED+"Source name as given in APEXHET fits files: "+" ".join(sourceList)+ENDCOLOR)
#                print(RED+"This probably means that the APEXHET fits files for the cube ARE NOT PRESENT,"+ENDCOLOR)
#                print(RED+"please download them from the archive and/or copy them to the current directory."+ENDCOLOR)
#                answ = str(input(BLUE+"In case the cube source name intentionally differs from the source name used in the APEXHET fits files, \n \
#                please indicate whether I should go on creating the Phase 3 files assuming all APEXHET files present belong to the cube? (y/n)"+ENDCOLOR))
#                answ = answ.lower()
#                if answ == 'y':
#                        print(GREEN+"Adding all APEXHET files despite discrepant source names!"+ENDCOLOR)
#                        for apexhet in fileList_tmp:
#                                with fits.open(apexhet) as hdul:
#                                        hetinf = hdul[0].header
#                                        fileList.append(apexhet)
#                                        extCont.append(int(hetinf['EXPTIME']))
#                                        proList.append(str(hetinf['HIERARCH ESO OBS PROG ID']))
#                else:
#                        print(RED+"Please make sure that all APEXHET fits files pertaining to the cube are present in the current directory!"+ENDCOLOR)
#                        sys.exit()
#
# else:
#        print()
#        print(RED+"Source name of cube and source names in APEXHET files do not (all) agree:"+ENDCOLOR)
#        print(RED+"Source name as given in cube: "+source+ENDCOLOR)
#        print(RED+"Source name(s) as given in APEXHET.fits files: "+" ".join(sourceList)+ENDCOLOR)
#        print()
#        if source.upper() in list(map(str.upper, sourceList)):
#                print(GREEN+"The source name as given in the cube is identical to the source name in a subset of the APEXHET fits files."+ENDCOLOR)
#                answ = str(input(BLUE+"I will "+RED+"ONLY"+BLUE+" include APEXHET files with source name identical to the cube. Is that ok? (y/n)"+ENDCOLOR))
#                answ = answ.lower()
#                if answ == 'y':
#                        for apexhet in fileList_tmp:
# with fits.open(apexhet) as hdul:
#                                        hetinf = hdul[0].header
#                                        sub_sou = str(hetinf['OBJECT'])
#                                        if source.upper() == sub_sou.upper():
#                                                fileList.append(apexhet)
#                                                extCont.append(int(hetinf['EXPTIME']))
#                                                proList.append(str(hetinf['HIERARCH ESO OBS PROG ID']))
#                else:
#                        answ = str(input(BLUE+"I will include "+RED+"ALL"+BLUE+" APEXHET files, assuming they all were used in a mosaic for the cube. Is that ok? (y/n)"+ENDCOLOR))
#                        answ = answ.lower()
#                        if answ == 'y':
#                                for apexhet in fileList_tmp:
#                                        with fits.open(apexhet) as hdul:
#                                                hetinf = hdul[0].header
#                                                fileList.append(apexhet)
#                                                extCont.append(int(hetinf['EXPTIME']))
#                                                proList.append(str(hetinf['HIERARCH ESO OBS PROG ID']))
#                        else:
#                                print(RED+"Please make sure that only APEXHET files pertaining to the present mosaic cube are present in the current directory"+ENDCOLOR)
#                                sys.exit()
#
#        else:
#                print(GREEN+"The source name as given in the cube doesn\'t match any of the source names in the APEXHET fits files."+ENDCOLOR)
#                answ = str(input(BLUE+'I will include '+RED+'ALL'+BLUE+' APEXHET files, assuming they all were used in a mosaic for the cube. Is that ok? (y/n) >>> '+ENDCOLOR))
#                answ = answ.lower()
#                if answ == 'y':
#                        for apexhet in fileList_tmp:
# with fits.open(apexhet) as hdul:
#                                        hetinf = hdul[0].header
#                                        fileList.append(apexhet)
#                                        extCont.append(int(hetinf['EXPTIME']))
#                                        proList.append(str(hetinf['HIERARCH ESO OBS PROG ID']))
#                else:
#                        print(RED+"Please make sure that only APEXHET files pertaining to the present mosaic cube are present in the current directory"+ENDCOLOR)
#                        sys.exit()
#
#### THE PREVIOUS SECTION WAS USED TO CHECK WHETHER OR NOT THE APEXHET FITS FILES PRESENT IN THE

fileList.sort()
print(
    GREEN
    + "These files will be associated with the cube:\n "
    + "\n".join(fileList)
    + ENDCOLOR
)

################################
# Number of files "APEXHET.*.fits" - TOTAL of files
################################
NrFilesBol = len(fileList)

if NrFilesBol >= 1:
    print("")
    print(GREEN + "Totat <APEXHET.*.fits> files >> " + ENDCOLOR, NrFilesBol)
    print(*fileList, sep="\t")
    print("")
else:
    print("")
    print(
        RED
        + "ESO ARCFILE data <APEXHET.*.fits> does not exist at the current Path"
        + ENDCOLOR
    )
    print(RED + "http://archive.eso.org/" + ENDCOLOR)
    sys.exit()

################################
# PROG_ID/PROGIDi more than one observing run - List modification.
################################
# print (proList)
proList = list(set(proList))  # set() doesn't have duplicate elements
lenProList = len(proList)

if lenProList == 1:
    prog_code = proList[0]
elif lenProList > 1:
    prog_code = "MULTI"
    # PROGIDi

print(prog_code)

################################
# Total integration time - EXPTIME
################################
"""
Total integration time per exposure is the sum of all individuals sub-integrations.

EXPTIME = TEXPTIME When:
The FULL map is not a mosaic. So the source verification could fails
because of the different in name and coverage, needs a different phase 3 treatment.
"""

ext_val = sum(extCont)

if ext_val < 5:
    print(
        RED
        + "There is an integration time inconsistency. Your ARCFILEs related to this program \
may be corrupted, Check your data"
        + ENDCOLOR
    )
    sys.exit()

################################
# Time notation for PROVi files
################################
for i in range(NrFilesBol):
    tmp_t = fileList[i]
    # tmp_t =  tmp_t.replace("_", ":")   # not needed, already in proper format in tap query
    # tmp_t =  tmp_t.replace(" ", ":")   # not needed, already in proper format in tap query
    fileTime = tmp_t[8:18] + " " + tmp_t[19:31]
    listTime.append(fileTime)
    # fileName = tmp_t  #[0:-5]   # not needed at all?
    # listName.append(fileName)   # not needed at all?
    print(tmp_t, fileTime)

################################
# Open the 1st file of ESO Science data <APEXHET.*.fits> for reference
################################
# NOT NEEDED, INFO WILL HAVE TO COME FROM TAP QUERY
# het_file = fileList[0]

# with fits.open(het_file) as hdul:
# het_info_0 = hdul[0].header
# het_info_1 = hdul[1].header

################################
# Extraction of Hierarch ESO info from data <APEXBOL.*.fits>.
################################
# ESO ID
# Scan technique
obs_tech = dprtech[0]
# obs_tech = het_info_0['HIERARCH ESO DPR TECH']
# obs_tech = 'SPECTRUM'
# Observation block ID, where available / Set of Observation Block (OB)
# ob_ID = het_info_0['HIERARCH ESO OBS ID'] # Not needed
# FEBE
# FEBE
##NEW for convenience I directly specify SUPERCAM-SCBE and SUPERCAM here, could also be taken from the get_frontBack() instead - or from the raw fits, ideally...
# febe = get_frontBack()
febe = "SUPERCAM-SCBE"  ##### THIS CAN'T BE TAKEN FROM THE TAP QUERY, IT SEEMS?
# febe	= het_info_0['FEBE1']
# instrument = 'SUPERCAM'
instrument = "APEXHET"
# febe	= 'HET230-XFFTS2'
# ind	= febe.index('-')
# front	= febe[0:ind]
# febe	= het_info_0['FEBE1']
# febe	= 'HET230-XFFTS2'
ind = febe.index("-")
front = febe[0:ind]
fac_yr = dicc[front][0]  # Year factor
factor = dicc[front][1]  # Conversion factor
fac_un = dicc[front][2]  # Uncertainty
# print (obs_tech, febe, front, factor)

################################
## Flux array
################################

print(
    CYAN
    + "APEX spectral line cubes may be ingested in the archive either on a \nbrightness temperature scale (K) or flux density (Jy) scale."
    + ENDCOLOR
)
print("")
print(CYAN + "In any case, a Jy/K conversion factor has to be provided" + ENDCOLOR)
print("")
print("")
print(
    CYAN
    + "- To obtain the conversion factor (Jy/K), please refer to the factor that best \
fits the time range of your observations."
    + ENDCOLOR
)
print("http://www.apex-telescope.org/telescope/efficiency/index.php")
print("")
print(
    CYAN
    + "Currently used Jy/K factor for "
    + str(front)
    + ": "
    + str(factor)
    + "+/-"
    + ""
    + str(fac_un)
    + "  (Year: "
    + str(fac_yr)
    + ")"
    + ENDCOLOR
)
print("")
answ = str(
    input(
        BLUE + "- Do you want to provide a different Jy/K factor? (y/n) >>> " + ENDCOLOR
    )
)
answ = answ.lower()

if answ == "y":
    try:
        factor = float(
            input(
                BLUE
                + "- Please, provide the new Jy/K factor (e.g. 40): >>>  "
                + ENDCOLOR
            )
        )
    except:
        print(RED + "Jy/K conversion factor must be a number" + ENDCOLOR)
        sys.exit()

    try:
        fac_un = float(
            input(BLUE + "- Please, provide the uncertainty : >>>  " + ENDCOLOR)
        )
    except:
        print(RED + "Jy/K uncertainty must be a number" + ENDCOLOR)
        sys.exit()

    print(GREEN + "Value set to: " + ENDCOLOR, factor, "+/-", fac_un)


try:
    if b_unit[0] == "K":
        print(
            CYAN
            + "- Your data product is currently expressed in terms of an antenna temperature scale [K]."
            + ENDCOLOR
        )
        print("")
        if b_unit != "K":
            print(
                CYAN
                + "The brightness temperature unit is given in "
                + b_unit
                + ", \npossibly no correction for beam efficiency has been applied?"
                + ENDCOLOR
            )
            print("")
            print(
                CYAN
                + "Please indicate whether a beam efficiency correction is still needed (if not, only the unit will be changed from '"
                + b_unit
                + "' to 'K')."
                + ENDCOLOR
            )
            answ = str(
                input(
                    BLUE
                    + " --- Do you want to apply a beam efficiency correction? (y/n) >>> "
                    + ENDCOLOR
                )
            )
            answ = answ.lower()
            if answ == "y":
                try:
                    beam_eff = float(
                        input(
                            BLUE
                            + "Please enter the beam efficiency (must be between 0 and 1 (greater than 0)) >>>"
                            + ENDCOLOR
                        )
                    )
                except:
                    print(
                        RED
                        + "Beam efficiency must be a number (between 0 and 1))"
                        + ENDCOLOR
                    )
                    sys.exit()

                if (beam_eff > 0.0) and (beam_eff <= 1.0):
                    print(
                        CYAN
                        + "        ... applying a beam efficiency of: "
                        + str(beam_eff)
                        + ENDCOLOR
                    )
                    print(
                        CYAN
                        + "        ... and changing brightness unit from '"
                        + b_unit
                        + "' to 'K'."
                        + ENDCOLOR
                    )
                    scidat = scidat / beam_eff
                    rmsdat = rmsdat / beam_eff
                else:
                    print(
                        RED
                        + "Beam efficiency must be between 0 and 1 (greater than 0 ...)"
                        + ENDCOLOR
                    )
                    sys.exit()

                print(beam_eff)

            else:
                print(
                    CYAN
                    + "        ... Changing brightness unit from '"
                    + b_unit
                    + "' to 'K' (without applying any further beam efficiency correction)"
                    + ENDCOLOR
                )

            b_unit = "K"
            b_unit_rms = "K"

        else:
            print(
                CYAN
                + "The brightness temperature unit is given in "
                + b_unit
                + ", probably a correction for beam efficiency has been applied?"
                + ENDCOLOR
            )
            answ = str(
                input(
                    BLUE
                    + "Please indicate whether this is correct (y) or whether you would like to apply an additional correction (n) >>> "
                    + ENDCOLOR
                )
            )
            answ = answ.lower()
            if answ != "y":
                try:
                    beam_eff = float(
                        input(
                            BLUE
                            + "Please enter the beam efficiency (must be between 0 and 1 (greater than 0)) >>>"
                            + ENDCOLOR
                        )
                    )
                except:
                    print(
                        RED
                        + "Beam efficiency must be a number (between 0 and 1))"
                        + ENDCOLOR
                    )
                    sys.exit()

                if (beam_eff > 0.0) and (beam_eff <= 1.0):
                    print(
                        CYAN
                        + "        ... applying a beam efficiency of: "
                        + str(beam_eff)
                        + ENDCOLOR
                    )
                    scidat = scidat / beam_eff
                    rmsdat = rmsdat / beam_eff
                else:
                    print(
                        RED
                        + "Beam efficiency must be between 0 and 1 (greater than 0 ...)"
                        + ENDCOLOR
                    )
                    sys.exit()

            else:
                print(
                    CYAN
                    + "        ... no further beam efficiency correction applied ..."
                    + ENDCOLOR
                )

        print("")
        print("")
        print(
            CYAN
            + "APEX spectral line cubes may be ingested in the archive either on a brightness temperature scale (K) or flux density (Jy) scale."
            + ENDCOLOR
        )
        answ = str(
            input(
                BLUE
                + "  --- Please indicate whether you want to convert the cube from 'K' to 'Jy' (y/n) >>> "
                + ENDCOLOR
            )
        )
        answ = answ.lower()

        if answ == "y":
            flux = scidat * factor
            flux_rms = rmsdat * factor
            b_unit = "Jy"
            b_unit_rms = "Jy"
        else:
            flux = scidat
            flux_rms = rmsdat

    elif b_unit == "mJy":
        print(
            CYAN
            + "- Your data product looks like it is calibrated in terms of Flux density [mJy]"
            + ENDCOLOR
        )
        answ = str(
            input(
                BLUE
                + "  --- Do you whish to convert to a brightness temperature (K) scale? (y/n) >>> "
                + ENDCOLOR
            )
        )
        answ = answ.lower()

        if answ == "y":
            flux = scidat / (1000.0 * factor)
            flux_rms = rmsdat / (1000.0 * factor)
            b_unit = "K"
            b_unit_rms = "K"
        else:
            flux = scidat
            flux_rms = rmsdat

    elif b_unit == "Jy":
        print(
            CYAN
            + "- Your data product looks like it is calibrated in terms of Flux density [Jy]"
            + ENDCOLOR
        )
        print("")
        answ = str(
            input(
                BLUE
                + "  --- Do you whish to convert to a brightness temperature (K) scale? (y/n) >>> "
                + ENDCOLOR
            )
        )
        answ = answ.lower()

        if answ == "y":
            flux = scidat / factor
            flux_rms = rmsdat / factor
            b_unit = "K"
            b_unit_rms = "K"
        else:
            flux = scidat
            flux_rms = rmsdat

    else:
        print(
            RED
            + "Your data seems to be neither in a brightness temperature (K) nor flux density (mJy/Jy) scale"
            + ENDCOLOR
        )
        print(
            RED
            + "Please check the calibration and provide a cube in proper flux/brightness units."
            + ENDCOLOR
        )
        sys.exit()

except:
    if b_unit[0] == "K":
        print(
            RED
            + "000 Intensity seems to be on a brightness temperature scale. Please try it once again."
            + ENDCOLOR
        )
        sys.exit()

    elif (b_unit == "Jy") or (b_unit == "mJy"):
        print(
            RED
            + "000 Intensity seems to be on a flux density scale. Please try it once again."
            + ENDCOLOR
        )
        sys.exit()

    else:
        print(
            RED
            + "000 Please check intensity units (should be on a 'K', 'Jy', or 'mJy' scale; currently it seems to be on a '"
            + b_unit
            + "' scale) and/or try again..."
            + ENDCOLOR
        )
        sys.exit()


def_val = np.nanmax(flux)
def_val = "%0.2f" % def_val
# print ('flux max', def_val)

## Determine median rms, convert to Jy if necessary (to fill BNOISE keyword)
rms_med = np.nanmedian(flux_rms)
if b_unit_rms[0] == "K":
    rms_med_Jy = rms_med * factor
else:
    rms_med_Jy = rms_med

## Collapse science cube to 2D white-light image
flux_white = np.nanmean(flux, axis=0)
# print(np.shape(flux_white))

## Frequency array from the WCS
## NEW the freqs array will be used to determine the minumum and maximum wavelength needed for the final header
freqs_tmp = freqArray(prihdr, specNr)
freqs = freqs_tmp[0]
refreq = freqs_tmp[1]
cdel = freqs_tmp[2]
# print(freqs, refreq, cdel)

# create HDUs holding the flux and rms cubes
hdu_fluxcube = fits.ImageHDU(data=flux, header=None, name="DATA_EXT")
hdu_rmscube = fits.ImageHDU(data=flux_rms, header=None, name="STAT_EXT")

## create HDU holding the white-light image
hdu_white = fits.PrimaryHDU(flux_white)

################################
# Validation - Time of observation
# Time in Julian day (start and stop)
################################
## NEW: This still needs to be worked out from the input fits files
time_1 = str(listTime[0])
time_2 = str(listTime[-1])

##Format example:
##Time_tmp1 = Time('2013-08-22 15:24:18', format='iso', scale='utc')
## Old APEXHET.2016-12-08T00_40_28.000.fits
## New APEXHET.2016-12-08T00:40:28.000
## rep APEXHET.2016-12-08T00_40_28.000

tmp_t1 = Time(time_1, format="iso", scale="utc")
tmp_t2 = Time(time_2, format="iso", scale="utc")

date_obs = time_1.replace(" ", "T")
time_star = float(round(tmp_t1.mjd, 5))
time_stop = float(round(tmp_t2.mjd, 5))
##Adding the last exposure time (s --> MJD)
time_stop = float(time_stop) + (extCont[-1] / (24 * 60 * 60))
time_stop = float("%0.5f" % time_stop) + 0.5  # --> seconds / + 0.5 uncertainty

##TEXPTIME (ext_val) cannot be greater than MJD-END - MJD-OBS (stop_start),
stop_start = (time_stop - time_star) * (24 * 60 * 60)  # --> seconds

# print(date_obs, time_star, time_stop, stop_start, ext_val)
if ext_val > stop_start:
    print("")
    print(RED + "Exposure time cannot be greater than (MJD-END - MJD-OBS)" + ENDCOLOR)
    print(
        str(ext_val) + "[s] ExpTime > " + str(int(stop_start)) + "[s] observing Time!"
    )
    sys.exit()


################################
### ESO - Programme validation .
################################
print("")
print(CYAN + "Entry section" + ENDCOLOR)
# prog_code

print(
    CYAN
    + ">>> "
    + ENDCOLOR
    + "Is this your ESO Programme Identification Code "
    + CYAN
    + str(prog_code)
    + ENDCOLOR
    + "?"
)
answ = str(input(" y / n >>>  "))
answ = answ.lower()


if answ == "y":
    print(GREEN + "ESO Programme Identification Code confirmed " + ENDCOLOR)

## NEW this probably needs to be adapted to include the new ESO Programme Code Format
elif answ == "n":
    proj_id = str(
        input("Enter the ESO Programme Identification Code / format:TP.C-NNNN(R) >>>  ")
    )
    chain_1 = "." in proj_id
    chain_2 = "-" in proj_id
    chain_3 = "(" in proj_id
    chain_4 = ")" in proj_id

    # print  (chain_1, chain_2, chain_3, chain_4)

    if chain_1 is True:
        if chain_2 is True:
            if chain_3 is True:
                if chain_4 is True:
                    prog_code = proj_id
                else:
                    print(
                        RED
                        + "The code  does not meet the minimum format requirements - EXIT"
                        + ENDCOLOR
                    )
                    sys.exit()
            else:
                print(
                    RED
                    + "The code  does not meet the minimum format requirements - EXIT"
                    + ENDCOLOR
                )
                sys.exit()
        else:
            print(
                RED
                + "The code  does not meet the minimum format requirements - EXIT"
                + ENDCOLOR
            )
            sys.exit()
    else:
        print(
            RED
            + "The code  does not meet the minimum format requirements - EXIT"
            + ENDCOLOR
        )
        sys.exit()

else:
    print(GREEN + "ESO Programme Identification Code keep it" + ENDCOLOR)

################################
# Bibliographic reference - Bibcode validation
################################
print("")
print(
    CYAN
    + ">>> "
    + ENDCOLOR
    + "Do you have Bibliographic reference for this data product?"
)
print(
    CYAN
    + "NOTE"
    + ENDCOLOR
    + ": REFERENC must be the 19-digit bibliograp otherwise press ENTER."
)
referen = str(input("Bibliographic reference: >>>  "))
ref_len = len(referen)  # (YYYYJJJJJVVVVMPPPPA)


try:
    if ref_len == 0:
        print(GREEN + "Bibcode entry left blank" + ENDCOLOR)
    elif ref_len > 0:
        ref_let = referen[-1]

        if ref_len == 19 and ref_let.isalpha():
            print(GREEN + "Bibcode keep it" + ENDCOLOR)
        else:
            print(RED + "That's invalid Bibcode, please try again." + ENDCOLOR)
            referen = str(input("Bibliographic reference (YYYYJJJJJVVVVMPPPPA): >>>  "))
            ref_len = len(referen)

            if ref_len == 0:
                print(GREEN + "Bibcode entry left blank" + ENDCOLOR)
            elif ref_len > 0:
                ref_let = referen[-1]
                if ref_len == 19 and ref_let.isalpha():
                    print(GREEN + "Bibcode keep it" + ENDCOLOR)
                else:
                    print(
                        RED
                        + "That's invalid Bibcode, Bibcode entry left blank"
                        + ENDCOLOR
                    )
                    referen = ""
except:
    pass

# ----------------------------------------------------------------------------------------
################################
# Keywords to fill
################################
crtDate = prihdr["DATE"]
crtoSofw = prihdr["ORIGIN"]
refrq4ap = refreq / 1.0e9  # GHz
# velo4xx  = velo  / 1000.0  # km/s

# Conversion of frequencies to wavelength (in vacuum)
WAVELMAX = (min(freqs) * u.GHz).to(u.nm, equivalencies=u.spectral()).value
WAVELMIN = (max(freqs) * u.GHz).to(u.nm, equivalencies=u.spectral()).value
WAVELMIN = float("%0.6f" % (WAVELMIN))
WAVELMAX = float("%0.6f" % (WAVELMAX))
# print(min(freqs), max(freqs))
# print(WAVELMIN, WAVELMAX)

## Average spectral resolving power
lambda_c = (WAVELMIN + WAVELMAX) / 2.0
bandwidth = WAVELMAX - WAVELMIN
channelwidth = bandwidth / (specNr - 1)
specres = lambda_c / channelwidth
# print(lambda_c, bandwidth, channelwidth, specNr)
# specres = float('%0.6f' %(specres))
print("Spectral resolving power: ", specres)

# SPEC_BIN keyword - not needed for cubes
# bin = (cdel * u.Hz).to(u.nm, equivalencies=u.spectral()).value
# bin  = float('%0.6f' %(bin))
##bin2 = WAVELMAX - WAVELMIN /(specNr - 1)
##print (cdel, bin, bin2)

##APERTURE keyword 7.8" * (800 / f [GHz])
## NEW this is actually not needed for cubes
D = 7.8 * (800 / refrq4ap)
D = ((D * u.arcsec).to(u.degree)).value
D = float("%0.8f" % (D))

## (effective) beam size, in arcsec
beam_size = 0.5 * (prihdr["BMAJ"] + prihdr["BMIN"]) * 3600.0
print(
    CYAN
    + "Effective beam size as obtained from the header of the input science cube is: "
    + str(beam_size)
    + ENDCOLOR
)
answ = str(
    input(
        BLUE
        + "   --- Would you like to enter a different beam size, e.g., to account for additional smoothing (y/n)? >>> "
        + ENDCOLOR
    )
)
answ = answ.lower()
if answ == "y":
    try:
        beam_size = float(input(BLUE + "Enter new beam size >>> " + ENDCOLOR))
    except:
        print(RED + "Beam size must be a number, please try again!" + ENDCOLOR)
        sys.exit()

print(prihdr)

# Delete keyword
del prihdr["BUNIT"]
del prihdr["DATAMIN"]  # NEW move to fluxcube extension instead?
del prihdr["DATAMAX"]  # NEW move to fluxcube extension instead?
del prihdr["LINE"]


equinox = 2000  # prihdr['EQUINOX']
c_frame = "FK5"  # het_info_0['RADECSYS']

if equinox == 2000 and c_frame == "FK5":
    coordRef = c_frame
else:
    print(
        RED
        + "Check your coord. frame. At APEX is FK5, why your it does not?"
        + ENDCOLOR
    )
    sys.exit()


# General
print("-------------------------Primary HDU-------------------------")
prihdr.set("PRODCATG", "SCIENCE.CUBE", "Data product category", 1)
prihdr.set("ORIGIN", "APEX", "Facility")
prihdr.set("TELESCOP", "APEX-12m", "Telescope name", after="ORIGIN")
# prihdr.set('INSTRUME', het_info_0['INSTRUME'], 'Instrument name' , after='TELESCOP')
prihdr.set("INSTRUME", instrument, "Instrument name", after="TELESCOP")
prihdr.set("FEBE1", febe, "APEX frontend/backend combination", after="INSTRUME")
prihdr.set("OBJECT", source, "Target designation", after="FEBE1")
prihdr.set("EXTEND", True, "Extensions are present")
# prihdr.set('RA'      , het_info_0['RA'], '[deg] Derived from native frame configuration', after='OBJECT')
# prihdr.set('DEC'     , het_info_0['DEC'],'[deg] Derived from native frame configuration', after='RA')
prihdr.set("RA", ra, "[deg] Derived from native frame configuration", after="OBJECT")
prihdr.set("DEC", dec, "[deg] Derived from native frame configuration", after="RA")
prihdr.set("EQUINOX", int(equinox), "Standard FK5 (years)", after="DEC")
prihdr.set("RADESYS", coordRef, "Coordinate reference frame", after="EQUINOX")
# prihdr.set('TIMESYS' , het_info_0['TIMESYS'], 'Time system for MJD', after='RADESYS')
prihdr.set("TIMESYS", "TAI", "Time system for MJD", after="RADESYS")
prihdr.set(
    "EXPTIME", ext_val, "[s] Total integration time per exposure", after="TIMESYS"
)
prihdr.set(
    "TEXPTIME", ext_val, "[s] Total integration time of all exposures", after="EXPTIME"
)
prihdr.set("MJD-OBS", time_star, "Start of observations (days)", after="TEXPTIME")
prihdr.set("MJD-END", time_stop, "End of observations (days)", after="MJD-OBS")
prihdr.set(
    "DATE-OBS", date_obs, "start of observation ISO 8601 format", after="MJD-END"
)
if prog_code == "MULTI":
    # prihdr.set('PROG_ID' , prog_code  , after='MJD-END') # NEW for the moment change to after TEXPTIME, needs to be changed back one MJD-END is given
    prihdr.set("PROG_ID", prog_code, after="TEXPTIME")
    for i in range(len(proList)):
        z = i + 1
        z = str(z)
        prihdr.set(
            "PROGID" + z, proList[i], "ESO programme identification", before="OBID1"
        )
else:
    # prihdr.set('PROG_ID' , prog_code  , 'ESO programme identification', after='MJD-END')
    prihdr.set(
        "PROG_ID", prog_code, "ESO programme identification", after="TEXPTIME"
    )  # NEW for the moment change to after TEXPTIME, needs to be changed back one MJD-END is given

prihdr.set(
    "SPEC_RES", specres, "Average spectral resolving power", after="PROG_ID"
)  # NEW after PROG_ID
prihdr.set("WAVELMIN", WAVELMIN, "[nm] Minimum wavelength", after="SPEC_RES")
prihdr.set("WAVELMAX", WAVELMAX, "[nm] Maximum wavelength", after="WAVELMIN")
## NEW not allowed for cubes: prihdr.set('SPEC_BIN', abs(bin),'[nm] Wavelength bin size', after='WAVELMAX')
prihdr.set("PROCSOFT", crtoSofw, "Reduction software/system", before="DATE")
prihdr.set("OBSTECH", obs_tech, "Technique of observation")
## NEW TBD: include check for MAPMODE (OTF/SPIRALRAS/SPIRAL)
prihdr.set("MAPMODE", "OTF", "APEX map mode")
prihdr.set("FLUXCAL", "ABSOLUTE", "Characterises the flux calibration", after="OBSTECH")
prihdr.set("BNOISE", rms_med_Jy, "Median rms (Jy)")

# List of Provenance files from lista
for i in range(NrFilesBol):
    b = i + 1
    b = str(b)
    prihdr.set(
        "PROV" + b, fileList[i].replace(" ", ":") + ".fits", "Original science file"
    )

prihdr.set("NCOMBINE", NrFilesBol, "Number of scans")
prihdr.set("JYFACTOR", factor, "[Jy/K] The Jansky to Kelvin conversion factor")
prihdr.set("JYUNIT", "[Jy/K]", "The Jansky to Kelvin conversion factor unit")
prihdr.set("JYUNCERT", fac_un, "[Jy/K] The Jansky to Kelvin conversion uncertainty")
prihdr.set("SKY_RES", beam_size, "effective beam size (arcsec)")
# prihdr.set('SKY_RERR', beam_size_error, 'error on beam size (arcsec)')
# prihdr.set('SKY_RERR', 0.0, 'error on beam size (arcsec)')
prihdr.set("REFERENC", referen, "Bibliographic reference")
# prihdr.set('ASSOCi'   , '' , 'Category of the associated non-FITS')
# prihdr.set('ASSONi'   , '' , 'List of files (FITS / non-FITS) associated to this data product.')

# print (prihdr)
# sys.exit()
################################
# Bintable
################################
## Creating the structure of the file: empty primary HDU followed by the extensions with the data and the rms.
phdu = fits.PrimaryHDU(header=prihdr)
## NEW change to use fluxcube instead of spectral table: hdulist = fits.HDUList([phdu, table_hdu])
hdulist = fits.HDUList([phdu, hdu_fluxcube, hdu_rmscube])

## Modify the extension header by adding the required Phase 3 keywords:
exthdr = hdulist[1].header
extrmshdr = hdulist[2].header

################################
# Populating science data extension header
################################
exthdr.set("HDUCLASS", "ESO", "class name (ESO format)")
exthdr.set("HDUDOC", "SDP", "ESO Science Data Products standard")
exthdr.set("HDUVERS", "SDP version 7", "version number")
exthdr.set("HDUCLAS1", "IMAGE", "data classification")
exthdr.set("HDUCLAS2", "DATA", "this extension contains the science data")
exthdr.set("ERRDATA", "STAT_EXT", "pointer to the error extension")

## propagating RA, DEC and OBJECT from the primary header. Such keywords need to be in the PHDU first.
exthdr["RA"] = (prihdr["RA"], "[deg] Spectroscopic target position (J2000.0)")
exthdr["DEC"] = (prihdr["DEC"], "[deg] Spectroscopic target position (J2000.0)")
exthdr["OBJECT"] = (prihdr["OBJECT"], "Target designation")

## propagating WCS from primary header to fluxcube extension

exthdr["NAXIS"] = prihdr["NAXIS"]
exthdr["NAXIS1"] = prihdr["NAXIS1"]
exthdr["NAXIS2"] = prihdr["NAXIS2"]
exthdr["NAXIS3"] = prihdr["NAXIS3"]
exthdr.set("BUNIT", b_unit, "", after="OBJECT")
exthdr["CTYPE1"] = prihdr["CTYPE1"]
exthdr["CRVAL1"] = prihdr["CRVAL1"]
## TEMP exthdr['CDELT1'] = prihdr['CDELT1']
exthdr["CRPIX1"] = prihdr["CRPIX1"]

## NEW check whether CUNITs are there and correctly set to 'deg'
if "CUNIT1" in prihdr:
    exthdr["CUNIT1"] = prihdr["CUNIT1"]
    if prihdr["CUNIT1"] == "":
        print(
            CYAN
            + "WCS first axis in input frame has no unit given (CUNIT1 is there, but empty), while a unit is required (and needs to be degrees)."
            + ENDCOLOR
        )
        if "CDELT1" in prihdr:
            pixelsize = prihdr["CDELT1"] * 3600.0
        else:
            if "CD1_1" in prihdr:
                pixelsize = prihdr["CD1_1"] * 3600.0
            else:
                print(
                    RED
                    + "... neither CDELT1 nor CD1_1 are in the input fits header. Please check the WCS headers of the input file"
                    + ENDCOLOR
                )
                sys.exit()

        print(
            CYAN
            + "Assuming that CDELT1 (or CD1_1) is in degrees, the pixel size is "
            + str(pixelsize)
            + " arcsec"
            + ENDCOLOR
        )
        answ = str(
            input(
                BLUE
                + "Unit of first WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> "
                + ENDCOLOR
            )
        )
        answ = answ.lower()
        if answ == "y":
            exthdr["CUNIT1"] = "deg"
        else:
            print(
                RED
                + "CUNITs for celestial coordinate system MUST be degrees."
                + ENDCOLOR
            )
            print(
                RED
                + "Please modify the WCS coordinate system of your input file(s) accordingly."
                + ENDCOLOR
            )
            sys.exit()
    else:
        if prihdr["CUNIT1"] != "deg":
            print(
                RED
                + "CUNITs for celestial coordinate system MUST be degrees."
                + ENDCOLOR
            )
            print(
                RED
                + "Please modify the WCS coordinate system of your input file(s) accordingly."
                + ENDCOLOR
            )
            sys.exit()

else:
    print(
        CYAN
        + "WCS first axis in input frame has no unit (CUNIT1) given, while a unit is required (and needs to be degrees)."
        + ENDCOLOR
    )
    if "CDELT1" in prihdr:
        pixelsize = prihdr["CDELT1"] * 3600.0
    else:
        if "CD1_1" in prihdr:
            pixelsize = prihdr["CD1_1"] * 3600.0
        else:
            print(
                RED
                + "... neither CDELT1 nor CD1_1 are in the input fits header. Please check the WCS headers of the input file"
                + ENDCOLOR
            )
            sys.exit()

    print(
        CYAN
        + "Assuming that CDELT1 (or CD1_1) is in degrees, the pixel size is "
        + str(pixelsize)
        + " arcsec"
        + ENDCOLOR
    )
    answ = str(
        input(
            BLUE
            + "Unit of first WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> "
            + ENDCOLOR
        )
    )
    answ = answ.lower()
    if answ == "y":
        exthdr.set("CUNIT1", "deg", "", after="CRPIX1")
    else:
        print(RED + "CUNIT for celestial coordinate system MUST be degrees." + ENDCOLOR)
        print(
            RED
            + "Please modify the WCS coordinate system of your input file(s) accordingly."
            + ENDCOLOR
        )
        sys.exit()

exthdr["CTYPE2"] = prihdr["CTYPE2"]
exthdr["CRVAL2"] = prihdr["CRVAL2"]
exthdr["CRPIX2"] = prihdr["CRPIX2"]

## NEW check whether CUNITs are there and correctly set to 'deg'
if "CUNIT2" in prihdr:
    exthdr["CUNIT2"] = prihdr["CUNIT2"]
    if prihdr["CUNIT2"] == "":
        print(
            CYAN
            + "WCS second axis in input frame has no unit given (CUNIT2 is there, but empty), while a unit is required (and needs to be degrees)."
            + ENDCOLOR
        )
        if "CDELT2" in prihdr:
            pixelsize = prihdr["CDELT2"] * 3600.0
        else:
            if "CD2_2" in prihdr:
                pixelsize = prihdr["CD2_2"] * 3600.0
            else:
                print(
                    RED
                    + "... neither CDELT2 nor CD2_2 are in the input fits header. Please check the WCS headers of the input file"
                    + ENDCOLOR
                )
                sys.exit()

        print(
            CYAN
            + "Assuming that CDELT2 (or CD2_2) is in degrees, the pixel size is "
            + str(pixelsize)
            + " arcsec"
            + ENDCOLOR
        )
        answ = str(
            input(
                BLUE
                + "Unit of second WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> "
                + ENDCOLOR
            )
        )
        answ = answ.lower()
        if answ == "y":
            exthdr["CUNIT2"] = "deg"
        else:
            print(
                RED
                + "CUNITs for celestial coordinate system MUST be degrees."
                + ENDCOLOR
            )
            print(
                RED
                + "Please modify the WCS coordinate system of your input file(s) accordingly."
                + ENDCOLOR
            )
            sys.exit()
    else:
        if prihdr["CUNIT2"] != "deg":
            print(
                RED
                + "CUNITs for celestial coordinate system MUST be degrees."
                + ENDCOLOR
            )
            print(
                RED
                + "Please modify the WCS coordinate system of your input file(s) accordingly."
                + ENDCOLOR
            )
            sys.exit()

else:
    print(
        CYAN
        + "WCS second axis in input frame has no unit (CUNIT2) given, while a unit is required (and needs to be degrees)."
        + ENDCOLOR
    )
    if "CDELT2" in prihdr:
        pixelsize = prihdr["CDELT2"] * 3600.0
    else:
        if "CD2_2" in prihdr:
            pixelsize = prihdr["CD2_2"] * 3600.0
        else:
            print(
                RED
                + "... neither CDELT2 nor CD2_2 are in the input fits header. Please check the WCS headers of the input file"
                + ENDCOLOR
            )
            sys.exit()

    print(
        CYAN
        + "Assuming that CDELT2 (or CD2_2) is in degrees, the pixel size is "
        + str(pixelsize)
        + " arcsec"
        + ENDCOLOR
    )
    answ = str(
        input(
            BLUE
            + "Unit of second WCS axis will be set to 'deg' (degree). Is that ok? (y/n) >>> "
            + ENDCOLOR
        )
    )
    answ = answ.lower()
    if answ == "y":
        exthdr.set("CUNIT2", "deg", "", after="CRPIX2")
    else:
        print(
            RED + "CUNITs for celestial coordinate system MUST be degrees." + ENDCOLOR
        )
        print(
            RED
            + "Please modify the WCS coordinate system of your input file(s) accordingly."
            + ENDCOLOR
        )
        sys.exit()

exthdr["CTYPE3"] = prihdr["CTYPE3"]
exthdr["CRVAL3"] = prihdr["CRVAL3"]
exthdr["CRPIX3"] = prihdr["CRPIX3"]

if "CD3_3" in prihdr:
    exthdr["CD3_3"] = prihdr["CD3_3"]
    if "CD1_3" in prihdr:
        exthdr["CD1_3"] = prihdr["CD1_3"]
    else:
        exthdr.set("CD1_3", 0.0, "Transformation matrix element", after="CD3_3")
    if "CD2_3" in prihdr:
        exthdr["CD2_3"] = prihdr["CD2_3"]
    else:
        exthdr.set("CD2_3", 0.0, "Transformation matrix element", after="CD1_3")
    if "CD3_1" in prihdr:
        exthdr["CD3_1"] = prihdr["CD3_1"]
    else:
        exthdr.set("CD3_1", 0.0, "Transformation matrix element", after="CD2_3")
    if "CD3_2" in prihdr:
        exthdr["CD3_2"] = prihdr["CD3_2"]
    else:
        exthdr.set("CD3_2", 0.0, "Transformation matrix element", after="CD3_1")
else:
    if "CDELT3" in prihdr:
        exthdr.set(
            "CD3_3",
            prihdr["CDELT3"],
            "Transformation matrix element spectral resolution",
            after="CRPIX3",
        )
        exthdr.set("CD1_3", 0.0, "Transformation matrix element", after="CD3_3")
        exthdr.set("CD2_3", 0.0, "Transformation matrix element", after="CD1_3")
        exthdr.set("CD3_1", 0.0, "Transformation matrix element", after="CD2_3")
        exthdr.set("CD3_2", 0.0, "Transformation matrix element", after="CD3_1")
    else:
        print(
            RED
            + "Can't find any information on spectral resolution (CDELT3 or CD3_3), this seems not to be a cube?"
            + ENDCOLOR
        )
        sys.exit()

if "CUNIT3" in prihdr:
    print(prihdr["CUNIT3"])
    exthdr["CUNIT3"] = prihdr["CUNIT3"]
else:
    print()
    print(RED + "CUNIT3 is not given." + ENDCOLOR)
    print(CYAN + "Spectral axis header values are:")
    print(CYAN + "Spectral axis type: " + prihdr["CTYPE3"] + ENDCOLOR)
    print(CYAN + "Reference channel value: " + str(prihdr["CRVAL3"]) + ENDCOLOR)
    print(CYAN + "Channel width: " + str(prihdr["CDELT3"]) + ENDCOLOR)
    if prihdr["CTYPE3"] == "VRAD":
        print(CYAN + "Frequency axis is in radial velocity units." + ENDCOLOR)
        answ = str(
            input(
                BLUE
                + "Unit of frequency axis will be set to 'm/s'. Is that ok? (y/n) >>> "
                + ENDCOLOR
            )
        )
        answ = answ.lower()
        if answ == "y":
            funits = "m/s"
        else:
            funits = str(
                input(
                    BLUE
                    + "Please enter units for the frequency axis (e.g., m/s, km/s) >>> "
                    + ENDCOLOR
                )
            )
            print("Frequency axis units will be set to: " + funits)
    else:
        funits = str(
            input(
                BLUE
                + "Please enter units for the frequency axis (e.g., m/s, km/s, GHz) >>> "
                + ENDCOLOR
            )
        )
        print("Frequency axis units will be set to: " + funits)

    exthdr.set("CUNIT3", funits, "", after="CRPIX3")


if (
    ("CD1_1" in prihdr)
    and ("CD1_2" in prihdr)
    and ("CD2_1" in prihdr)
    and ("CD2_2" in prihdr)
):
    print("CDi_j transformation matrix is present")
else:
    print("CDi_j transformation matrix is not present in input header.")
    if "CROTA2" in prihdr:
        if ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
            print(
                "Computing CDi_j transformation matrix elements from CDELTi and CROTA2 keywords"
            )
            cd12 = -1.0 * prihdr["CDELT2"] * np.sin(prihdr["CROTA2"])
            cd21 = prihdr["CDELT1"] * np.sin(prihdr["CROTA2"])
            # print(cd12, cd21)
            if cd12 == -1.0 * cd12:
                cd12 = np.abs(cd12)
            if cd21 == -1.0 * cd21:
                cd21 = np.abs(cd21)
            # print(cd12, cd21)
            exthdr.set(
                "CD1_1",
                prihdr["CDELT1"] * np.cos(prihdr["CROTA2"]),
                "Transformation matrix element",
            )
            exthdr.set("CD1_2", cd12, "Transformation matrix element")
            exthdr.set("CD2_1", cd21, "Transformation matrix element")
            exthdr.set(
                "CD2_2",
                prihdr["CDELT2"] * np.cos(prihdr["CROTA2"]),
                "Transformation matrix element",
            )
        else:
            print(
                RED
                + "Neither transformation matrix CDi_j nor CDELT1/CDELT2 are given"
                + ENDCOLOR
            )
            print(
                RED + "Please check the WCS system of your input data cube." + ENDCOLOR
            )
            sys.exit()
    else:
        if ("CDELT1" in prihdr) and ("CDELT2" in prihdr):
            print("Neither transformation matrix nor CROTA2 keyword found.")
            print(
                "Creating transformation matrix assuming that there is no rotation to the coordinate system"
            )
            exthdr.set("CD1_1", prihdr["CDELT1"], "Transformation matrix element")
            exthdr.set("CD1_2", 0.0, "Transformation matrix element")
            exthdr.set("CD2_1", 0.0, "Transformation matrix element")
            exthdr.set("CD2_2", prihdr["CDELT2"], "Transformation matrix element")
        else:
            print(
                RED
                + "Neither transformation matrix CDi_j nor CDELT1/CDELT2 are given"
                + ENDCOLOR
            )
            print(
                RED + "Please check the WCS system of your input data cube." + ENDCOLOR
            )
            sys.exit()

exthdr["SPECSYS"] = prihdr["SPECSYS"]
exthdr["RESTFREQ"] = prihdr["RESTFREQ"]
exthdr["VELO-LSR"] = prihdr["VELO-LSR"]
exthdr["IMAGFREQ"] = prihdr["IMAGFREQ"]

################################
# Populating rms/error extension header
################################

extrmshdr.set("HDUCLASS", "ESO", "class name (ESO format)")
extrmshdr.set("HDUDOC", "SDP", "ESO Science Data Products standard")
extrmshdr.set("HDUVERS", "SDP version 7", "version number")
extrmshdr.set("HDUCLAS1", "IMAGE", "data classification")
extrmshdr.set("HDUCLAS2", "ERROR", "this extension contains the error data")
extrmshdr.set("HDUCLAS3", "RMSE", "error type")
extrmshdr.set("SCIDATA", "DATA_EXT", "pointer to the science extension")

## propagating RA, DEC and OBJECT from the primary header. Such keywords need to be in the PHDU first.
extrmshdr["RA"] = (prihdr["RA"], "[deg] Spectroscopic target position (J2000.0)")
extrmshdr["DEC"] = (prihdr["DEC"], "[deg] Spectroscopic target position (J2000.0)")
extrmshdr["OBJECT"] = (prihdr["OBJECT"], "Target designation")

## propagate WCS from primary header to rms cube extension
## original fits cubes have additional CROTAn keywords which are not allowed, should be converted into a transformation matrix instead (but are 0 for my case)
extrmshdr["NAXIS"] = prihdr["NAXIS"]
extrmshdr["NAXIS1"] = prihdr["NAXIS1"]
extrmshdr["NAXIS2"] = prihdr["NAXIS2"]
extrmshdr["NAXIS3"] = prihdr["NAXIS3"]
extrmshdr.set("BUNIT", b_unit_rms, "", after="OBJECT")
extrmshdr["CTYPE1"] = prihdr["CTYPE1"]
extrmshdr["CRVAL1"] = prihdr["CRVAL1"]
## TEMP extrmshdr['CDELT1'] = prihdr['CDELT1']
extrmshdr["CRPIX1"] = prihdr["CRPIX1"]
extrmshdr["CUNIT1"] = exthdr["CUNIT1"]
extrmshdr["CTYPE2"] = prihdr["CTYPE2"]
extrmshdr["CRVAL2"] = prihdr["CRVAL2"]
## TEMP extrmshdr['CDELT2'] = prihdr['CDELT2']
extrmshdr["CRPIX2"] = prihdr["CRPIX2"]
extrmshdr["CUNIT2"] = exthdr["CUNIT2"]
extrmshdr["CTYPE3"] = prihdr["CTYPE3"]
extrmshdr["CRVAL3"] = prihdr["CRVAL3"]
extrmshdr["CRPIX3"] = prihdr["CRPIX3"]
# extrmshdr['CDELT3'] = prihdr['CDELT3'] # use CD3_3 instead
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

################################
# Populating white-light image header
################################
whitehdr = hdu_white.header

whitehdr.set("PRODCATG", "ANCILLARY.IMAGE.WHITELIGHT")
whitehdr.set("NAXIS", 2, "")
## propagating RA, DEC and OBJECT from the primary header. Such keywords need to be in the PHDU first.
whitehdr["RA"] = (prihdr["RA"], "[deg] Spectroscopic target position (J2000.0)")
whitehdr["DEC"] = (prihdr["DEC"], "[deg] Spectroscopic target position (J2000.0)")
whitehdr["OBJECT"] = (prihdr["OBJECT"], "Target designation")

## propagate WCS from primary header to fluxcube extension

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
####################################

# Delete some primary headers that are not needed any more (or were moved to the etension header)
del prihdr["CTYPE*..."]
del prihdr["CRPIX*..."]
del prihdr["CRVAL*..."]
del prihdr["CDELT*..."]
del prihdr["CROTA*..."]
del prihdr["CUNIT1*..."]
del prihdr["NAXIS"]
del prihdr["NAXIS1"]
del prihdr["NAXIS2"]
del prihdr["NAXIS3"]
del prihdr["RESTFREQ"]
del prihdr["VELO-LSR"]
del prihdr["IMAGFREQ"]

print(prihdr)
print("")
print(exthdr)
print("")
print(extrmshdr)
print("")
print(whitehdr)
# ----------------------------------------------------------------------------------------

################################
# "Virtual" output files names
################################
outFluname = source + "_3DcubePh3.fits"
outWhitename = source + "_3DcubePh3_whitelight.fits"
prihdr.set("ASSON1", outWhitename)
if attlmv == "y":
    outLMVname = source + "_3Dcube_GILDASformat.lmv"
    prihdr.set("ASSON2", outLMVname)
    prihdr.set("ASSOC2", "ANCILLARY.CUBE.LMV")
    md5 = hashlib.md5(open(lmvIn, "rb").read()).hexdigest()
    prihdr.set("ASSOM2", md5)

# ----------------------------------------------------------------------------------------
print(
    GREEN
    + "New FITS files with your calibrated data product are created!. \n \
They contain the data (flux and rms cube) in a HDU structure and with a header compliant with the Phase 3 requirements.\
In addition a white-light image and a fits file containing the flux cube with the original HDU structure is created. \n \
The data should then be ready to be included in a Phase 3 submission."
    + ENDCOLOR
)
if attlmv == "y":
    print(
        GREEN
        + "An IRAM GILDAS format .lmv cube has been associated and will also have to be provided in the Phase 3 submission."
        + ENDCOLOR
    )
print(GREEN + "Fingers crossed..." + ENDCOLOR)
# ----------------------------------------------------------------------------------------
################################
# Files Creation
################################
hdulist[0].header = prihdr
hdulist[1].header = exthdr
hdulist[2].header = extrmshdr

hdu_white.header = whitehdr

###Flux - Intensity - Signal
# Create file out.fits containing an HDU constructed from data and
# header containing both CHECKSUM and DATASUM cards.
# fits.writeto( outFluname, scidat, prihdr, overwrite=True, checksum=True) #with checksun & datasum &scidata
# print (nameFig)
# ----------------------------------------------------------------------------------------

## Saving the cube and rms into a new file
hdulist.writeto(outFluname, checksum=True, overwrite=True)
print(outFluname + " created...")

## Saving white light image into a new file
hdu_white.writeto(outWhitename, checksum=True, overwrite=True)
print(outWhitename + " created...")

## Creating copy of GILDAS format lmv cube for upload:
if attlmv == "y":
    os.system("cp -i " + lmvIn + " " + outLMVname)
    print(outLMVname + " created (copy of original GILDAS .lmv format cube)")
