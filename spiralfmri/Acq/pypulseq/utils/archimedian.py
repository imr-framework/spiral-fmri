
import numpy as np
import math

def archimedian_spiral(smax, gmax, T, N, FOV, rmax):
    """Creates a variable density spiral waveform with trajectory
    k(t) = r(t) exp(i*q(t))
    Where q is the same as theta, r and q are chosen to satisfy:
    1) Maximum gradient amplitudes and slew rates.
    2) Maximum gradient due to FOV, where FOV can vary with k-space radius r, as
    FOV(r) = F0 + F1*r + F2*r*r

    Parameters
    ----------
    smax : float
        Maximum slew rate in mT/m
    gmax : float
        Maximum gradient it T/m/s
    T : float
        Gradient raster time in seconds
    N : int
        Number of spiral shots
    FOV :
        FOV coefficient with respect to r in meters
    rmax : float
        Maximum k-space radius in 1/m

    Returns
    -------
    k : numpy.ndarray
        K-space trajectory (kx + iky) in 1/m
    g : numpy.ndarray
        Gradient waveform (Gx + iGy) in Hz/m

    """

    if N < 1:
        raise ValueError('Number of shots has to be at least 1')
    if rmax == 0:
        raise ValueError('Maximum k-space radius cannot be 0')
    elif rmax < 0:
        rmax = - rmax
    if FOV <= 0:
        raise ValueError('FOV cannot be zero or negative')

    F0 = FOV
    F1 = 0
    F2 = 0

    gamma = 42576000 # Hz/T
    gmax = gmax * gamma * 1e-3 # Hz/m
    smax = smax * gamma # Hz/m/s

    oversamp = 8
    To = T / oversamp


    q0 = 0
    q1 = 0
    theta = np.zeros((1000000,))
    r = np.zeros((1000000,))
    r0 = 0
    r1 = 0

    time = np.zeros((1000000,))
    t = 0
    count = 0

    while r0 < rmax:
        q2, r2 = q2r21(smax, gmax, r0, r1, To, T, N, [F0, F1, F2])

        q1 = q1 + q2 * To
        q0 = q0 + q1 * To
        t = t + To

        r1 = r1 + r2 * To
        r0 = r0 + r1 * To

        count += 1
        theta[count] = q0
        r[count] = r0
        time[count] = t

    r = r[np.arange(int(oversamp / 2)-1, count + 1, oversamp)]
    # r = r[np.arange(int(oversamp / 2), count + oversamp, oversamp)]
    theta = theta[np.arange(int(oversamp / 2) - 1, count + 1, oversamp)]
    time = time[np.arange(int(oversamp / 2) - 1, count + 1, oversamp)]

    ltheta = 4 * math.floor(len(theta) / 4)
    r = r[:ltheta]
    # r = r[:ltheta + 1]
    theta = theta[:ltheta]
    time = time[:ltheta]

    k = r * np.exp(1j * theta)

    g = (np.hstack((k, 0)) - np.hstack((0, k))) / T
    # 1 / gamma * (np.hstack((k, 0)) - np.hstack((0, k))) / T
    g = g[:len(k)]

    s = (np.hstack((g,0)) - np.hstack((0,g))) / T
    s = s[:len(k)]

    return k, g #, s, time, r, theta

def q2r21(smax, gmax, r, r1, T, Ts, N, Fcoeff):

    gamma = 42576000  # Hz/T
    F = 0
    dFdr = 0
    for rind in range(1, len(Fcoeff) + 1):
        F = F + Fcoeff[rind -1] * r ** (rind - 1)
        if rind > 1:
            dFdr = dFdr + (rind - 1) * Fcoeff[rind - 1] * r ** (rind - 2)

    GmaxFOV = 1 / F / Ts # Hz/m
    # GmaxFOV = 1 / gamma / F / Ts
    Gmax = min(GmaxFOV, gmax)

    maxr1 = math.sqrt((Gmax) ** 2 / (1 + (2 * math.pi * F * r / N) ** 2)) # Hz/m
    # math.sqrt((gamma * Gmax) ** 2 / (1 + (2 * math.pi * F * r / N) ** 2))

    if r1 > maxr1:
        r2 = (maxr1 - r1) / T
    else:
        twopiFoN = 2 * math.pi * F / N
        twopiFoN2 = twopiFoN ** 2

        A = 1 + twopiFoN2 * r * r
        B = 2 * twopiFoN2 * r * r1 * r1 + 2 * twopiFoN2 / F * dFdr * r * r * r1 * r1
        C = twopiFoN2 ** 2 * r * r * r1 ** 4 + 4 * twopiFoN2 * r1 ** 4 + (2 * math.pi / N * dFdr) ** 2 * r * r * r1 ** 4 + 4 * twopiFoN2 / F * dFdr * r * r1 ** 4 - smax ** 2

        rts = qdf(A, B, C)
        r2 = rts[0].real

        slew = 1 / (r2 - twopiFoN2 * r * r1 ** 2 + 1j * twopiFoN * (2 * r1 ** 2 + r * r2 + dFdr / F * r * r1 ** 2))
        # 1 / gamma * (r2 - twopiFoN2 * r * r1 ** 2 + 1j * twopiFoN * (2 * r1 ** 2 + r * r2 + dFdr / F * r * r1 ** 2)) # Hz/m/s

        sr = abs(slew) / smax

        if sr > 1.01:
            raise ValueError(f'Slew rate violation {sr * 100}')

    q2 = 2 * math.pi / N * dFdr * r1 ** 2 + 2 * math.pi * F / N * r2

    return q2, r2

def qdf(a, b, c):

    d = b ** 2 - 4 * a * c
    r1 = (-b + math.sqrt(d)) / (2 * a)
    r2 = (-b - math.sqrt(d)) / (2 * a)

    return r1, r2