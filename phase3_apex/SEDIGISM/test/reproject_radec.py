from astropy.io import fits
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp

hdu = fits.open('G301_HC3N_Tmb_DR1.fits')[0]
hdu_2d = hdu
hdu_2d.data = hdu.data[0, :, :]
hdu_2d.header['NAXIS'] = 2
del hdu_2d.header['*3*']

wcs_out, shape_out = find_optimal_celestial_wcs(hdu_2d, frame='fk5', auto_rotate=True)
hdr_out = wcs_out.to_header()

array, footprint = reproject_interp(hdu_2d, wcs_out.to_header(), shape_out=shape_out)