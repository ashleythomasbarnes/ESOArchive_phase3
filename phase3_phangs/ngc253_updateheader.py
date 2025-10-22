from astropy.io import fits 

hdu_cube = fits.open("NGC0253.fits")
hdu_hdr  = fits.open("NGC0253_header.fits")

hdu_cube[0].header = hdu_hdr[0].header
hdu_cube[1].header = hdu_hdr[1].header
hdu_cube[2].header = hdu_hdr[2].header
hdu_cube[3].header = hdu_hdr[3].header

hdu_cube.writeto("NGC0253.fits", overwrite=True, checksum=True)

hdu_cube.close()
hdu_hdr.close()

hdu_cube_check = fits.open("NGC0253.fits")

print(hdu_cube_check[0].header)
print(hdu_cube_check[1].header)
print(hdu_cube_check[2].header)
print(hdu_cube_check[3].header)