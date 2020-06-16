"""
Author: Marina Manso Jimeno
Last Updated: 07/11/2019
"""

import numpy as np

def g2k(gx, gy, gz):
    """
    Transforms from gradient to k-space trajectory

    Parameters
    ----------
    gx : numpy.ndarray
        Gradient waveform along the x-axis
    gy : numpy.ndarray
        Gradient waveform along the y-axis
    gz : numpy.ndarray
        Gradient waveform along the z-axis

    Returns
    -------
    kx : numpy.ndarray
        k-space trajectory along the x-axis
    ky : numpy.ndarray
        k-space trajectory along the x-axis
    kz : numpy.ndarray
        k-space trajectory along the x-axis
    """

    g = []
    for i in range(len(gx)):

        g.append(complex(gx[i],gy[i]))
    gts = 10e-6 # gradient sample duration

    ktemp = np.cumsum(g)*gts # cycles/m
    kx = np.real(ktemp)
    ky = np.imag(ktemp)
    kz = np.cumsum(gz) * gts  # cycles/m

    return kx, ky, kz