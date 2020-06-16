'''
Helper functions for spiral-fmri reconstruction
Author : Marina Manso Jimeno
Last updated : 06/11/2020
'''

import numpy as np
import math

def data_check_file_inputs(raw_dat: np.ndarray, ktraj: np.ndarray):
    '''
    Checks that the dimensions of the files inputted for raw data and k-space trajectory match

    Parameters
    ----------
    raw_dat : numpy.ndarray
         k-space raw data from the scanner with dimensions [Npoints x kz x Nframes x Nchannels]
    ktraj : numpy.ndarray
        k-space trajectory with dimensions [Npoints x Nframes x kz x 3]
    '''

    if raw_dat.shape[0] != ktraj.shape[0] or raw_dat.shape[1] != ktraj.shape[1] or raw_dat.shape[2] != ktraj.shape[2]:
        raise ValueError('The raw data and k-space trajectory provided do not match')

def data_check_getparams_input(raw_dat: np.ndarray, params: dict):
    '''
    Checks that the dimensions of the file inputted for raw data and the parameters from getparams() match

    Parameters
    ----------
    raw_dat : numpy.ndarray
         k-space raw data from the scanner with dimensions [Npoints x kz x Nframes x Nchannels]
    params : dict
        Acquisition parameters dictionary from getparams() function
    '''
    kz = params['acq_params']['nz']
    nframes = params['fmri']['nframes']
    if kz != raw_dat.shape[1]:
        raise ValueError('Acquisiton parameters and raw data file do not match: number of stacked spirals (kz)')
    if nframes != raw_dat.shape[2]:
        raise ValueError('Acquisiton parameters and raw data file do not match: number of frames')

def rfphs_compensate(raw_dat : np.ndarray, rfphs_angle : float):
    '''
    Compensates for the RF phase offset applied during acquisition

    Parameters
    ----------
    raw_dat : numpy.ndarray
        k-space raw data from the scanner with dimensions [Npoints x kz x Nframes x Nchannels]
    rfphs_angle : float
        RF phase angles in degrees

    Returns
    -------
    rfphs_corr_dat : numpy.ndarray
        Compensated data
    '''
    rfphs = 0
    rf_spoil_seed_cnt = 0
    rf_spoil_seed = rfphs_angle

    kz = raw_dat.shape[1]
    nframes = raw_dat.shape[2]
    total_spirals = kz * nframes

    rfphs_v = np.zeros((total_spirals, 1))
    for i in range(len(rfphs_v)):
        rfphs_v[i] = rfphs
        rfphs = rfphs + (rf_spoil_seed / 180 * math.pi) * rf_spoil_seed_cnt
        rf_spoil_seed_cnt += 1

    rfphs_v = np.roll(rfphs_v, 1) # check if it's negative or positive

    count = 0
    rfphs_corr_kspace = np.zeros(raw_dat.shape,  dtype=complex)
    for fr in range(nframes):
        for z in range(kz):
            rfphs_corr_kspace[:, z, fr, :] = raw_dat[:, z, fr, :] * np.exp(-1j * rfphs_v[count])
            count += 1

    return rfphs_corr_kspace

def spiral_alignment(raw_dat : np.ndarray, nleafs: int):
    '''
    Re-order the raw data to correctly align the spiral shots from different frames for sliding window reconstruction

    Parameters
    ----------
    raw_dat : numpy.ndarray
        k-space raw data from the scanner with dimensions [Npoints x kz x Nframes x Nchannels]
    nleafs : int
        Number of spiral leafs of a fully-sampled slice

    Returns
    -------
    aligned_dat : numpy.ndarray
        k-space raw data arranged combining leafs from nleafs consecutive time points [Npoints x nleafs x kx x Nframes - nleafs + 1 x Nchannels]
    '''
    kz = raw_dat.shape[1]
    nframes = raw_dat.shape[2]

    rot_fr = np.arange(1,nleafs + 1)

    aligned_dat = np.zeros((raw_dat.shape[0], nleafs, kz, nframes - nleafs + 1, raw_dat.shape[-1]), dtype=complex)
    for fs_ksp in range(aligned_dat.shape[3]):
        rot_sl = rot_fr
        for sl in range(kz):
            for leaf in range(nleafs):
                aligned_dat[:, leaf, sl, fs_ksp, :] = raw_dat[:, sl, fs_ksp + np.where(rot_sl == leaf + 1)[0][0], :]
                '''aligned_dat[:, 1, sl, fs_ksp, :] = raw_dat[:, sl, fs_ksp + np.where(rot_sl == 2)[0][0], :]
                aligned_dat[:, 2, sl, fs_ksp, :] = raw_dat[:, sl, fs_ksp + np.where(rot_sl == 3)[0][0], :]'''
            rot_sl = np.roll(rot_sl, -1)
        rot_fr = np.roll(rot_fr, -1)

    return aligned_dat