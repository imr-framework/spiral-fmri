from math import sqrt, ceil

import numpy as np


def trapwave2(area, mxg, mxs, rasterTime):
    '''
    Create gradient trapezoid with given area.

    Parameters
    ----------
    area : float
        Area in mT/m*sec
    mxg : int
        Maximum gradient limit in mT/m
    mxs : int
        Maximum slew rate limit  in T/m/s
    rasterTime : float
        Gradient raster time in seconds

    Returns
    -------
    waveform : numpy.ndarray
        Trapezoidal gradient waveform in mT/m
    '''
    # area = area*1e3
    # Do all calculations as positive then flip at end if negative

    if rasterTime <= 0:
        raise ValueError('Raster time has to be larger than 0')

    if area < 0:
        area = -area
        invertarea = 1
    elif area == 0:
        raise ValueError('Area cannot be 0')
    else:
        invertarea = 0

    # reduce peak amp/slew so it passes hardware checks
    mxg = 0.995 * mxg
    mxs = 0.995 * mxs * 1e3 # mT/m/s

    # construct trapezoid
    dt = rasterTime # seconds
    tr = mxg/mxs # seconds # time for ramp to full scale
    Acrit = mxs * tr ** 2 # mT/m*s
    dg = mxs * dt # mT / m
    if area < Acrit:
        rtime = sqrt(area / mxs) # seconds
        n = ceil(rtime / dt)
        ramp = np.arange(n + 1).reshape(1, n+1) * dg
        waveform = np.concatenate((ramp, np.fliplr(ramp)), axis=1)
    else:
        nr = ceil(tr / dt)
        ramp = np.arange(nr) * mxs * dt
        areaRamps = 2 * np.sum(ramp) *dt
        np1 = ceil((area-areaRamps) / mxg / dt)
        plat = mxg * np.ones((1,np1))
        waveform = np.concatenate((ramp.reshape(1,len(ramp)), plat, hflip(ramp)), axis=1)

    # scale down to desired area
    wavArea = np.sum(waveform) * dt
    if wavArea < area:
        raise ValueError('Can"t scale down to desired area. bug in code')

    waveform = waveform / wavArea * area

    return waveform

def hflip(arr, cp = 0):
    '''
    Flips the ramp of a trapezoidal gradient

    Parameters
    ----------
    arr : numpy.ndarray
        array

    Returns
    ------
    newarr : numpy.ndarray
        flipped array
    '''

   # s = arr.shape
    s2 = arr.shape[0]

    if cp == 0:
        cp = (s2 + 1) / 2

    cp = cp + s2
    temp = np.concatenate((arr, arr, arr), axis=0)
    temp = temp.reshape(1, len(temp))

    newarr = 0 * arr.reshape(1, s2)

    for k in range(s2):
        newarr[0, k] = temp[0, int(2 * cp - k - s2) - 2]

    return newarr
