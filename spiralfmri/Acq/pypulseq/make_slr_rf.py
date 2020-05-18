# Copyright of the Board of Trustees of Columbia University in the City of New York

"""
Radio-Frequency SLR pulse generation
Author: Marina Manso Jimeno
Last Updated: 04/20/2020
"""

from math import sqrt, log10, log, atan2, cos, sin, pi

import numpy as np
from scipy.signal import firls
from scipy.signal import resample
from scipy.signal import resample_poly

from pypulseq.opts import Opts

from spiralfmri.Acq.pypulseq.makeSystemlength import makeSystemlength
from spiralfmri.Acq.pypulseq.trapwave2 import trapwave2

def make_slr_rf(flip_angle: float, slice_thickness: float, time_bw_product: float, duration: float, ncycles: int,
                    system: dict):
    '''
    Creates a slice-selective SLR pulse with gradient crusher (or balancing blip) before it.

    Parameters
    ----------
    flip_angle : float
        Flip angle in degrees.
    slice_thickness : float
        Slab thickness in millimeters (mm) of accompanying slice select trapezoidal event. The slice thickness determines the area of the
        slice select event.
    time_bw_product : float
        Time-bandwidth product of SLR pulse
    duration : float
        Duration in seconds (s).
    ncycles : int
        number of cycles of phase (spoiler) (0 = balanced)
    system : dict
        System limits. Default is a system limits object initialised to default values.

    Returns
    -------
    rf : SimpleNamespace
        Radio-frequency slr pulse event.
    gz : SimpleNamespace, optional
        Accompanying slice select trapezoidal gradient event. Returned only if `slice_thickness` is not zero.
    '''
    if flip_angle < 0 or flip_angle > 360:
        raise ValueError('Flip angle in degrees has to be within [0, 360]')
    if slice_thickness <= 0:
        raise ValueError('Slice thickness has to be larger than 0')
    if time_bw_product < 1.5:
        raise ValueError('Time bandwidth product has to be at least 1.5')
    if duration <= 0:
        raise ValueError('Pulse duration cannot be 0 or negative')

    # TODO: Remove commented lines
    # TODO: Remove frequency offset from input
    pulse_type = 'ex' # Only supported type by now
    # make RF waveform
    dt = system['rf_raster_time']# s
    res = round(duration / dt) # Number of (10us) samples in RF waveform
    duration = res * dt

    nrf = 200 # number of samples for SLR design
    rf = JP_dzrf(nrf, time_bw_product, pulse_type)
    rf = np.real(rf)
    rf = resample_poly(rf, res, nrf)


    flip = pi/2

    rf = flip * rf / np.sum(rf, axis=0)
    rf = rf / dt / 2 / pi # Hz

    rf = np.concatenate((rf.reshape(len(rf),1), np.zeros((len(rf) % 2, 1))), axis=0) # Make duration even
    npix = len(rf)

    iref,_ = np.where(rf == np.max(rf)) # Center of RF pulse

    # make slice-select gradient waveform
    bw = time_bw_product / duration # Hz
    gplateau = bw / (system['gamma'] * slice_thickness) * 1e3 # mT / m

    if  gplateau > system['max_grad']:
        raise ValueError('Gradient exceeds the maximum gradient')

    # slice-select trapezoid
    gss = gplateau * np.ones((1,npix))
    s = system['max_slew'] * 1e3 * dt * 0.999 # mT/m
    gss_ramp = np.arange(s, gplateau, s)
    if np.all(gss_ramp == 0):
        gss_ramp = np.asarray([[0]])

    if gplateau - gss_ramp[-1] > 0:
        gss_ramp = np.concatenate((gss_ramp.reshape(1, len(gss_ramp)), ((gplateau + gss_ramp[-1])/2).reshape(1,1)), axis=1)

    gss_trap = np.concatenate((gss_ramp, gss, np.fliplr(gss_ramp)), axis=1)
    iref = iref[0] + gss_ramp.shape[1]

    # slice-select rephaser trapezoid

    arearep = np.sum(gss_trap[:,iref:]) * dt  # mT / m * s
    gzrep = -trapwave2(arearep, system['max_grad'], system['max_slew'], dt)

    # slice-select prephaser trapezoid

    areaprep = np.sum(gss_trap[:,:iref]) * dt # mT / m * s
    gzprep = -trapwave2(areaprep, system['max_grad'], system['max_slew'], dt)


    #irep = (np.concatenate((gzprep, gss_trap), axis=1)).shape[1]
    #iref = iref + gzprep.shape[1]
    gex = np.concatenate((gzprep, gss_trap, gzrep), axis=1) * 1e-3 * system['gamma'] # Hz/m
    #idep = gzprep.shape[1]

    # make gss and rf the same length
    rf = np.vstack((0 * gzprep.reshape(-1,1), np.zeros((gss_ramp.shape[1],1)), rf, np.zeros((gss_ramp.shape[1] + gzrep.shape[1],1))))
    # rf = np.concatenate((0 * gzprep, np.zeros((1, gss_ramp.shape[1])), rf.T, np.zeros((1, gss_ramp.shape[1] + gzrep.shape[1]))), axis=1)

    # ensure that duration is on a 40 us (4 sample) boundary
    rf = makeSystemlength(rf, system['grad_raster_time'])
    gex = makeSystemlength(gex, system['grad_raster_time'])

    # make sure waveforms start and end at zero
    rf = np.vstack((0, rf, 0))
    # rf = np.concatenate((np.asarray([[0]]), rf, np.asarray([[0]])), axis=1)
    gex = np.vstack((0, gex.reshape(-1,1), 0))
    # gex = np.concatenate((np.asarray([[0]]), gex, np.asarray([[0]])), axis=1)

    rf = rf / 90 * flip_angle

    # slice offset frequency
    # gfreq = system['gamma']*gplateau*freq_offset



    #remove -0 values
    # rf[np.where(rf == -0.)] = 0
    # gex[np.where(gex == -0.)] = 0
    return rf, gex



def JP_dzrf(np, tb, ptype: str = 'ex', ftype: str = 'ls', d1 = 0, d2 = 0):
    '''
    Designs an rf pulse

    Parameters
    ----------
    np : integer
        Number of points
    tb:
        Time-bandwidth product
    ptype : string, optional
        Pulse type (ex = pi/2 excitation pulse)
    ftype : string, optional
        Filter design method (ls = least squares)
    d1 : integer, optional
        Passband ripple. Default is 0.
    d2 : integer, optional
        Stopband ripple. Default is 0.

    Returns
    -------
    rf : numpy.ndarray
        rf waveform
    '''
    # J. Pauly dzrf.m
    d1 = 0.01
    d2 = d1

    bsf = sqrt(1 / 2)
    d1 = sqrt(d1 / 2)
    d2 = d2 / sqrt(2)

    b = JP_dzls(np, tb, d1, d2)
    b = bsf * b
    rf = JP_b2rf(b)

    return rf

def JP_dzls(nf,tb,d1,d2):
    '''
    Designs a least squares filter.

    Parameters
    ----------
    nf : int
        Filter length
    tb : float
        Time-bandwidth
    d1 : float
        Pass band ripple
    d2 : float
        Stop band ripple

    Returns
    -------
    h : numpy.ndarray
        Least squares filter
    '''
    di = JP_dinf(d1,d2)
    w = di / tb
    f = np.asarray([0, (1 - w) * (tb / 2), (1 + w) * (tb / 2), nf / 2]) / (nf / 2)
    m = np.asarray([1, 1, 0, 0])
    w = np.asarray([1, d1 / d2])

    h = firls(nf-1, f, m, w)

    h = resample(h, nf)

    return h

def JP_dinf(d1,d2):
    '''
    Calculates D infinity for a linear phase filter
    
    Parameters
    ----------
    d1 : float
        Pass band ripple
    d2 : float
        Stop band ripple

    Returns
    -------
    d : float
        D infinity for a linear phase filter
    '''
    
    a1 = 5.309e-3
    a2 = 7.114e-2
    a3 = -4.761e-1
    a4 = -2.66e-3
    a5 = -5.941e-1
    a6 = -4.278e-1

    l10d1 = log10(d1)
    l10d2 = log10(d2)

    # Right now only dealing with floats, size = 1,1
    # m1, n1 = l10d1.shape
    m1 = 1
    n1 = 1
    if  m1 < n1:
        l10d1 = l10d1.T

    # m2, n2 = l10d2.shape
    m2 = m1
    n2 = n1
    if m2 < n2:
        l10d2 = l10d2.T

    # l = len(d2)

    # Cannot transpose a float
    # d = (a1 * l10d1 * l10d1 + a2 *l10d1 +a3) * l10d2.T + (a4 * l10d1 * l10d1 + a5 * l10d1 + a6) * np.ones((1,1))
    d = (a1 * l10d1 * l10d1 + a2 * l10d1 + a3) * l10d2 + (a4 * l10d1 * l10d1 + a5 * l10d1 + a6) #* np.ones(1, 1)
    return d

def JP_b2rf(bc):
    '''
    This function takes a beta polynomial and returns an RF pulse
    A minimum phase alpha polynomial is intermediately computed,
    followed by the SLR transform.

    Works for complex RF pulses

    Parameters
    ----------
    bc : numpy.ndarray
        beta polynomial coefficients

    Returns
    -------
    rf : numpy.ndarray
        RF pulse waveform
    '''

    ac = JP_b2a(bc)
    rf = JP_ab2rf(ac, bc)

    return rf

def JP_b2a(bc):
    '''
    This function takes a b polynomial, and returns the minimum phase,
    minimum power a polynomial

    Parameters
    ----------
    bc : numpy.ndarray
        Beta polynomial coefficients

    Returns
    -------
    aca : numpy.ndarray
        Minimum phase alpha polynomial
    '''

    n = len(bc)
    bcp = bc
    bl = len(bc)
    blp = bl * 8
    bcp_add = np.zeros((blp-len(bcp)))
    bcp = np.concatenate((bcp, bcp_add), axis=0)
    bf = np.fft.fft(bcp)
    bfmax = np.max(np.abs(bf))
    if bfmax >= 1:
        bf = bf / (1e-8 + bfmax)
    afa = JP_mag2mp(np.sqrt(1 - bf * np.conj(bf)))
    aca = np.fft.fft(afa) / blp
    aca = aca[np.arange(n-1,-1,-1)]
    aca = np.reshape(aca, bc.shape)

    return aca

def JP_mag2mp(x):
    '''
    Take the magnitude of the fft of a signal, and return the
    fft of the analytic signal.

    Parameters
    ----------
    x : numpy.ndarray
        magnitude of analytic signal fft

    Returns
    -------
    a : numpy.ndarray
        fft of analytic signal
    '''

    n = len(x)
    xl = np.log(x)
    xlf = np.fft.fft(xl)
    xlfp = np.zeros(xlf.shape, dtype=complex)
    xlfp[0] = xlf[0]
    xlfp[1:int(n/2)] = 2 * xlf[1:int(n/2)]
    xlfp[int(n/2)] = xlf[int(n/2)]
    xlfp[int(n/2)+1:] = 0 * xlf[int(n/2)+1:]
    xlaf = np.fft.ifft(xlfp)
    a = np.exp(xlaf)

    return a

def JP_ab2rf(ac,bc):
    '''
    Takes two polynomials for alpha and beta, and return an
    rf waveform that would generate them.

    Parameters
    ----------
    ac : numpy.ndarray
        polynomial for alpha
    bc : numpy.ndarray
        polynomial for beta

    Returns
    -------
    rf2 : numpy.ndarray
        rf waveform that produces alpha and beta under the hard pulse approximation
    '''
    n = len(ac)
    j = 1j
    c = np.zeros((n,))
    s = np.zeros((n,),dtype=complex)
    rf = np.zeros((n,),dtype=complex)
    for i in np.arange(n-1, -1, -1):
        c[i] = sqrt(1 / (1 + abs(bc[i] / ac[i]) ** 2))
        s[i] = np.conj(c[i] * bc[i] / ac[i])
        theta = atan2(abs(s[i]), c[i])
        psi = np.angle(s[i])
        rf[i] = 2 * (theta * cos(psi) + j * theta * sin(psi))
        acn = c[i] * ac + s[i] * bc
        bcn = -np.conj(s[i]) * ac + c[i] * bc

        ac = acn[1:i+1]
        bc = bcn[:i]


    return rf