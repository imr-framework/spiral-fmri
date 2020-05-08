import numpy as np

from Acq.trapwave2 import trapwave2

def make_balanced(g, max_grad, max_slew, system):
    '''
    Adds a trapezoid at end of the gradient waveform to make the total area zero.

    Parameters
    ----------
    g : numpy.ndarray
        Gradient waveform in Hz/m
    system : dict
        System limits and specifications

    Returns
    -------
    gbal : numpy.ndarray
        Balanced gradient in Hz/m
    '''
    if len(g.shape) == 1:
        g = g.reshape(g.shape[0],1)
    elif g.shape[1] > 1 and g.shape[0] == 1:
        g = g.reshape(-1,1)
    elif g.shape[0] > 1 and g.shape[1] > 1:
        raise ValueError('Only one dimensional waveforms are allowed')

    if all([ v == 0 for v in g ]):#g.all() == 0:
        raise ValueError('The gradient has to have nonzero elements')

    if 'gamma' not in system or 'grad_raster_time' not in system:
        raise ValueError('Gamma and gradient raster time have to be specified')

    max_slew_Hzms = max_slew * system['gamma']  # Hz/m/s

    # dt  seconds
    dt = system['grad_raster_time']
    # ramp to zero
    dg = - np.sign(g[-1]) * max_slew_Hzms * 0.995 * dt # Hz/sample
    ramp = np.arange(g[-1], 0+dg, dg)
    g = np.concatenate((g, ramp[:-1].reshape(len(ramp) - 1, 1)))

    # change of units
    g = g * 1e3 / system['gamma'] # mT/m

    # add balancing trapezoid
    area = np.sum(g) * dt # mT * s / m
    gblip = trapwave2(np.abs(area), max_grad * 0.995, max_slew * 0.995, dt) * 1e-3 * system['gamma'] # Hz/m
    g = g * 1e-3 * system['gamma']  # Hz/m
    gbal = np.concatenate((g, -np.sign(area) * gblip.T))

    return gbal