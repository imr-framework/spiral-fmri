"""
Functions for RF and Gradient raster times interpolation

For RF:
GE raster time = 4 micro seconds
Siemens raster time = 1 micro second

For Gradients:
GE raster time = 4 micro seconds
Siemens raster time = 10 micro seconds

Author: Marina Manso Jimeno
Last Updated: 07/11/2019
"""
import numpy as np

def rf_interpolate(rf_signal, rf_raster_time):
    GE_rf_raster_time = 4e-6
    T = len(rf_signal) * GE_rf_raster_time  # pulse duration
    t_GE = np.linspace(0, T - GE_rf_raster_time, T / GE_rf_raster_time)
    t_out = np.linspace(0, T - rf_raster_time, T / rf_raster_time)
    rf_out = np.interp(t_out, t_GE, rf_signal.flatten())

    return rf_out


def grad_interpolate(grad_signal, grad_raster_time):
    GE_grad_raster_time = 4e-6
    T = len(grad_signal) * GE_grad_raster_time  # pulse duration
    t_GE = np.linspace(0, T - GE_grad_raster_time, T / GE_grad_raster_time)
    t_out = np.linspace(0, T - grad_raster_time, T / grad_raster_time)
    grad_out = np.interp(t_out, t_GE, grad_signal.flatten())

    return grad_out