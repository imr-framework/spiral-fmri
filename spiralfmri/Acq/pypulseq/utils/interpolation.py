"""
Author: Marina Manso Jimeno
Last Updated: 05/06/2020
"""
import numpy as np

from spiralfmri.Acq.pypulseq.utils.makeSystemlength import makeSystemlength

def rf_interpolate(rf_signal, rf_raster_time=1e-6):
    '''
    RF waveform interpolation to meet raster time requirements

    Parameters
    ----------
    rf_signal : numpy.ndarray
        RF signal with raster time 4 microseconds
    rf_raster_time : float
        Raster time for interpolation. Defalult is 1 microsecond

    Returns
    -------
    rf_out :  numpy.ndarray
        Interpolated RF waveform with new raster time
    '''
    if rf_raster_time <= 0:
        raise ValueError('Raster time has to be larger than 0')

    GE_rf_raster_time = 4e-6
    T = len(rf_signal) * GE_rf_raster_time  # pulse duration
    t_GE = np.linspace(0, T - GE_rf_raster_time, T / GE_rf_raster_time)
    t_out = np.linspace(0, T - rf_raster_time, T / rf_raster_time)
    rf_out = np.interp(t_out, t_GE, rf_signal.flatten())
    rf_out = makeSystemlength(rf_out, rf_raster_time)

    return rf_out


def grad_interpolate(grad_signal, grad_raster_time=10e-6):
    '''
    Gradient waveform interpolation to meet raster time requirements

    Parameters
    ----------
    grad_signal : numpy.ndarray
        Gradient signal with raster time 4 microseconds
    grad_raster_time : float
        Raster time for interpolation. Default is 10 microseconds

    Returns
    -------
    grad_out : numpy.ndarray
        Interpolated gradient waveform with new raster time
    '''
    if grad_raster_time <= 0:
        raise ValueError('Raster time has to be larger than 0')

    GE_grad_raster_time = 4e-6
    T = len(grad_signal) * GE_grad_raster_time  # pulse duration
    t_GE = np.linspace(0, T - GE_grad_raster_time, T / GE_grad_raster_time)
    t_out = np.linspace(0, T - grad_raster_time, T / grad_raster_time)
    grad_out = np.interp(t_out, t_GE, grad_signal.flatten())
    grad_out = makeSystemlength(grad_out, grad_raster_time)

    return grad_out