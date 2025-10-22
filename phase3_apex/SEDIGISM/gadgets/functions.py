#!/usr/bin/env python
# Version: 2022-05-20

import numpy as np

##########################################
#BLACK	 = "\033[30m"
RED      = "\033[31m"		# For errors
GREEN    = "\033[32;3m"		# For outputs
#YELLOW   = "\033[33;1m"	# Not in use
#BLUE	 = "\033[34;1m"		# Not in use
MAGENTA  = "\033[35m"		# Header
CYAN     = "\033[36m"		# Main  words
#WHITE    = "\033[37m"
ENDCOLOR = "\033[0m"
##########################################

def header():
	print(MAGENTA+'##########################################'+ENDCOLOR)
	print (MAGENTA+'#APEX-Phase 3 Tool'+ENDCOLOR)
	print (MAGENTA+'#Creation of 1D FITS spectra in'+ENDCOLOR)
	print (MAGENTA+'#compliant with ESO Phase 3 requirements'+ENDCOLOR)
	print (MAGENTA+'#ESO-APEX 12m telescope'+ENDCOLOR)
	print (MAGENTA+'#version 3.2 / date 2022-Jun-20'+ENDCOLOR)	
	print (MAGENTA+'#contact: pvenegas@eso.org'+ENDCOLOR)
	#print (MAGENTA+'#       : xxx@apex-telescope.org'+ENDCOLOR)
	print (MAGENTA+'#Author : P. Venegas'+ENDCOLOR)
	print(MAGENTA+'##########################################'+ENDCOLOR)
	print('')
	return()

def info():
	##Modified!##### Check this info before realease
	print(RED+'Information:'+ENDCOLOR, '\n',\
	'* Only the RAW Science Data files (one per scan number) together with your Data Product will be used to produce a New Header,',\
	'calibration data will not be included.','\n',\
	'* The pertinent information of the RAW science files going to be query from the TAP database. This will creates an ascii',\
	'file (extension .tap) that going to be used as input of this Tool.','\n',\
	'* Therefore, only your Data Product is needed in this Folder to get started.')
	return()

	
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
 
	target=source

	print()
	print("Querying the ESO TAP service at %s" %(ESO_TAP_OBS))

	# --------------------------------------------------
	# The actual position of the selected target
	# is queried by the from_name() function,
	# which queries the CDS SESAME service
	# (http://cdsweb.u-strasbg.fr/cgi-bin/Sesame).
	# --------------------------------------------------

	query="""SELECT dp_id, exposure, prog_id, object, dp_tech, instrument, ra, dec 
	from dbo.raw
	where dp_id like 'APEXHET.%%'
  		and object = '%s'
  		and dp_cat='SCIENCE'""" % (target);

	print()
	print("Query:")
	print(query)

	res = tapobs.search(query=query,maxrec=1000)

	print()
	print(res.to_table())
	print()


	print()
	print("A total of "+str(len(res.to_table()))+" records were found matching the provided criteria.")

	table = res.to_table()

	filename = "%s.tap" % target
	print()
	print ("Results has been written to: %s" % filename)
	
	with open(filename, 'w') as f:
	 print (tabulate(table), file=f)
     
	f.close()
	return(table)


def instrList(myList):
	import sys
	c=0
	for i in range(len(myList)): #print((*myList), sep = "["+ +"]\n") / enumerate
		c = c+1
		print (CYAN+ myList[c-1]," ["+str(c)+"]" +ENDCOLOR)
	try:	
		febe_tmp = int(input(' >>> '))
	except :
		print (RED+'Backend is need it, try once again'+ENDCOLOR )
		sys.exit()
	
	febe_usr = myList[febe_tmp-1]	
	print( ' > '+ febe_usr)
	print ('')
	return(febe_usr)


def instrument(op_inst):
	print ('> Choose the option of your used Instrument-Backend [?]')
	if op_inst == 1 : #SEPIA	
		myList = ['GARD180-PBE_F','GARD180-XFFTS','NOVA660-PBE_F','NOVA660-XFFTS','SEPIA180-FFTS1','SEPIA180-PBE_F','SEPIA180-XFFTS','SEPIA345-FFTS1','SEPIA660-FFTS1','SEPIA660-PBE_F','SEPIA660-XFFTS']
		febe_usr = instrList(myList) 

	if op_inst == 2 : #NFLASH
		myList = ['NFLASH230-FFTS1','NFLASH460-FFTS1']
		febe_usr = instrList(myList) 

	if op_inst == 3 : #FLASH
		myList = ['FLASH345-AFFTS','FLASH345-FFTS1','FLASH345-FFTS4G','FLASH345-PBE_A','FLASH345-PBE_F','FLASH345-XFFTS','FLASH460-AFFTS','FLASH460-FFTS','FLASH460-FFTS1','FLASH460-PBE_A','FLASH460H-XFFTS','FLASH460L-FFTS1','FLASH460L-PBE_F','FLASH460L-XFFTS','FLASH810-FFTS','FLASH810-FFTS1','FLASH810-FFTS2','FLASH810-PBE_A']
		febe_usr = instrList(myList) 

	if op_inst == 4 : #PI230
		myList = ['PI230-FFTS4G','PI230-PBE_C','PI230-PBE_E']
		febe_usr = instrList(myList) 

	if op_inst == 5 : #SHFI
		myList = ['HET230-PBE_A','HET230-FFTS1','HET230-XFFTS2','HET345-PBE_A','HET345-FFTS1','HET345-FFTS2','HET345-XFFTS2','HET460-PBE_A','HET460-FFTS1','HET460-XFFTS2','HET1300-PBE_A','HET1300-FFTS1']
		febe_usr = instrList(myList) 	
		
	return(febe_usr)
	
		
def get_frontBack():	
	import sys
	print (CYAN+'Instrument / FEBE, i.e. HET230-XFFTS2'+ENDCOLOR)
	#print ('http://archive.eso.org/wdb/wdb/eso/apex/query')
	#print ('> Instrument used')
	print ('> This implementation supports the following heterodyne instruments:')
	print ('> Choose an option: '+CYAN+'SEPIA [1], nFLASH [2], FLASH [3], PI230 [4] and SHeFI [5]'+ENDCOLOR)
	op_inst = input('> Integer option >>> ')
	print("")
	OPTIONS = (1, 2, 3, 4, 5)
	
	try :
		op_inst = int(op_inst)
	except :
		print (RED+'Frontend is need it, try once again'+ENDCOLOR )
		sys.exit()
		
	if op_inst in OPTIONS : 
		febe_usr = instrument(op_inst)
	else:
		print (RED+'FEBE need it, try once again'+ENDCOLOR )
		sys.exit()
							
	return (febe_usr)
    

#Integer 
def is_integer_num(n):
    if isinstance(n, int):
        return True
    if isinstance(n, float):
        return n.is_integer()
    return False
    
    
#Plotting function.    
def plotFig1(source, front, wav, flux, b_unit):
	################################
	import matplotlib.pyplot as plt	
	################################
	fig = plt.figure()
	ax = fig.add_subplot(111)
	fig.subplots_adjust(top=0.9)
	fig.suptitle(str(front)+  '/ Source: ' + str(source), fontsize=12, fontweight='bold')
	ax.set_title('(Plot based on Jy-to-K factor provided)', fontsize=8)
	plt.title('This is a verification plot, not need it by Ph3 process', loc='right', y=-0.01 ,fontsize=6.5, fontstyle='italic', color='red')

	plt.plot(wav, flux)
	plt.xlabel('Frequency (GHz)')
	plt.ylabel('Intensity ('+str(b_unit)+')')
	plt.grid(True)
	nameFig = (str(source)+'_plotPh3.png')
	plt.savefig(nameFig, format = 'png')
	plt.show()
	plt.close()

	return(nameFig)  
	

def error_prop(fac_un, factor, noise, peak ):
	#Error propagation.
	#Flux_err key
	if np.isnan(peak) == False:
		E1 = (fac_un/factor)*100
		E2 = (noise/peak)*100
		fluxerr = ( (E1)**2  + (E2)**2 ) ** 0.5
		
	else:
		#E1 = (fac_un/factor)*100
		#fluxerr = ( (E1)**2  + (noise)**2 ) ** 0.5	#Noise + relative factor error.
		#fluxerr = ( noise * factor)	#If noise unit is K
		fluxerr = ( noise )				#if noise unit is Jy
	
	return (fluxerr)
	

def waveArray(prihdr, specNr):	 
	## Wavelength array from the WCS
	w0 		= prihdr['CRVAL1']
	cdel 	= prihdr['CDELT1']
	crpix1 	= prihdr['CRPIX1']
	refreq	= prihdr['RESTFREQ']
	velo	= prihdr['VELO-LSR']

	p0 = crpix1 * cdel + refreq
	wav = []
	print ('Frecuency axis in construction')
	for i in range(specNr):
  		w1 =  (p0 - (i) * cdel) / 1.0e9   # == 10**9, GHz
  		#print (w1)
  		wav.append(w1)
	return (wav, refreq, cdel)	