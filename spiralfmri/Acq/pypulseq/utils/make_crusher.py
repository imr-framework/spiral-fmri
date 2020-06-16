from spiralfmri.Acq.pypulseq.utils.trapwave2 import trapwave2
from spiralfmri.Acq.pypulseq.utils.makeSystemlength import makeSystemlength

import numpy as np

def make_crusher(ncycles, opslthick, gzarea, max_grad, max_slew, sys_lims):
    '''
    Creates a crusher gradient within the system limits given the number of phase cycles across a thickness

    Parameters
    ----------
    ncycles : int
        Number of cycles of phase across slice/slab
    opslthick : float
        Slice/slab thickness in meters
    gzarea : int
        Half-area of slice-select gradient in mT/m*sec. Default: 0.
    max_grad : int
        Max gradient in mT/m
    max_slew : float
        Max slew rate  in T/m/s
    sys_lims : dict
        System limits

    Returns
    -------
    gcrush : numpy.ndarray
        Crusher gradient waveform in Hz/m
    '''
    if ncycles <= 0:
        raise ValueError('Number of phase cycles cannot be equal or smaller than 0.')
    if opslthick == 0:
        raise ValueError('Thickness cannot be 0')
    if 'gamma' not in sys_lims:
        raise ValueError('Please specify gamma in the system limits')

    gamma = sys_lims['gamma'] # Hz/T

    area = ncycles / (gamma * opslthick) * 1e3 # mT/m*s

    dt = sys_lims['grad_raster_time'] # seconds
    gcrush = trapwave2(area - gzarea, max_grad, max_slew, dt) * 1e-3 * gamma # Hz/m
    gcrush = makeSystemlength(gcrush.reshape(-1,1), sys_lims['grad_raster_time'])
    gcrush[np.where(gcrush == -0.)] = 0
    return gcrush