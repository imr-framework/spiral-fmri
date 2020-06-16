'''
Reconstruction script for stack of spirals trajectories
Author: Marina Manso Jimeno
Last updated: 06/04/2020
'''

import numpy as np
import numpy.fft as fft

from spiralfmri.Recon.Python.utils.NUFFT_recon import nufft_init

def recon_SoS(raw_dat : np.ndarray, ktraj : np.ndarray, N : int):
    '''
    Parameters
    ----------
    raw_dat : numpy.ndarray
        Raw data array with shape (Npoints, Nleafs, Nz, Nchannels)
    ktraj : numpy.ndarray
        K-space trajectory coordinates array with shape (Npoints, Nleafs)
    N : int
        Image size

    Returns
    -------
    recon_vol_sos : numpy.ndarray
        Reconstructed image volume after channel combination using sum of squares (N, N, Nz)
    '''

    # Shape checks
    if raw_dat.shape[0] != ktraj.shape[0]:
        raise ValueError('Number of points from raw data and trajectory do not match')
    if raw_dat.shape[1] != ktraj.shape[1]:
        raise ValueError('Number of spiral leafs from raw data and trajectory do not match')

    if len(raw_dat.shape) < 4:
        raw_dat = np.expand_dims(raw_dat, axis=2)

    nufft_args = (N, 2 * N, 6) # Nd, Kd, Jd
    nufft_obj = nufft_init(ktraj, nufft_args)

    # Do IFT along kz
    if raw_dat.shape[2] > 1:
        raw_dat = fft.fftshift(fft.ifft(fft.ifftshift(raw_dat, axes=2), axis=2), axes=2)

    Nslices = raw_dat.shape[2]
    Nchannels = raw_dat.shape[-1]
    recon_vol = np.zeros((N, N, Nslices, Nchannels), dtype=complex)
    for z in range(Nslices):
        for ch in range(Nchannels):
            recon_vol[:, :, z, ch] = nufft_obj.solve(np.squeeze(raw_dat[:, :, z, ch]).flatten(), solver='cg', maxiter=50)

    # Sum of squares
    recon_vol_sos = np.sqrt(np.sum(abs(recon_vol) ** 2, 3))
    recon_vol_sos = np.fliplr(np.rot90(recon_vol_sos, -1))
    return recon_vol_sos