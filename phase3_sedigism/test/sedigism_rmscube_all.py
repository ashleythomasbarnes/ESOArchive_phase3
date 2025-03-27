import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
from tqdm import tqdm
import glob


def listfiles(fname,path,keep_path=False):

    pfilenames = glob.glob('%s%s' % (path, fname))
    filenames = []
    for pfilename in pfilenames:
        if keep_path:
            filenames.append(pfilename)
        else:
            filenames.append(pfilename[len(path)::])			
        
    return filenames


path = './'

fnames = listfiles('./data/*_DR1.fits', path, keep_path=True)

for fname in tqdm(fnames):

    # print('%s/%s' % (i+1len(fnames)))

    hdu = fits.open(fname)[0]

    hdu.data = hdu.data.squeeze()
    sizes = hdu.data.shape
    
    nchans = 100 # Number of line-free channels to consider
    free1 = hdu.data[sizes[0]-(nchans+1):sizes[0]-1,:,:]
    free2 = hdu.data[0:nchans,:,:]

    free = np.concatenate([free2,free1],axis=0)

    # rmsmap = np.nanstd(free, axis = 0)
    rmsmap = mad_std(free, axis=0, ignore_nan=True)

    mask = free > rmsmap *2
    free_masked = free.copy()
    free_masked[mask] = np.nan

    rmsmap = mad_std(free_masked, axis=0, ignore_nan=True)

    rmscube = np.repeat(rmsmap[np.newaxis,:,:], hdu.data.shape[0], axis=0)

    rmscube = fits.PrimaryHDU(rmscube, hdu.header)
    rmscube.writeto(fname.split('.fits')[0]+'_RMSCUBE.fits', overwrite=True)