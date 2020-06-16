'''
Helper methods for Non-Uniform Fast Fourier transform reconstruction
Author: Marina Manso Jimeno
Last updated: 06/04/2020
'''
import math
import numpy as np
import pynufft

def nufft_init(kt : np.ndarray, params):
    '''Initializes the Non-uniform FFT object

    Parameters
    ----------
    kt : numpy.ndarray
        K-space trajectory
    params : list
        Sequence parameters.

    Returns
    -------
    NufftObj : pynufft.linalg.nufft_cpu.NUFFT_cpu
        Non-uniform FFT Object for non-cartesian transformation
    '''
    kt_sc = math.pi / abs(np.max(kt))
    kt = kt * kt_sc

    om = np.zeros((kt.shape[0] * kt.shape[1],2))
    om[:, 0] = np.real(kt).flatten()
    om[:, 1] = np.imag(kt).flatten()

    NufftObj = pynufft.NUFFT_cpu()  # Create a pynufft object
    Nd = (params[0], params[0])  # image size
    Kd = (params[1], params[1])  # k-space size
    Jd = (params[2], params[2])  # interpolation size

    NufftObj.plan(om, Nd, Kd, Jd)
    return NufftObj