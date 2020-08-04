
import numpy as np

def getparams():
    """
    Acquisition parameters for 3D stack-of-spirals fMRI scan

    Returns
    -------
    params : dict
        Dictionary containing the acquisition and specific sequence parameters
    """

    # system limits
    gamma = 42576000 # Hz/T
    sysSiemens = {'gamma': gamma, 'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                 'grad_raster_time': 10e-6, 'rf_raster_time': 1e-6}
    #TODO: Try having maximum grad as 32 for GE
    sysGE = {'gamma': gamma, 'max_grad': 31, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s', 'grad_raster_time': 4e-6, 'rf_raster_time': 4e-6 }

    # Do not modify system parameters for now. These are tested and work.
    if sysSiemens['max_grad'] != 32 or sysSiemens['grad_unit'] != 'mT/m':
        raise ValueError('Please do not modify Siemens gradient maximum value or units.')
    if sysSiemens['max_slew'] != 130 or sysSiemens['slew_unit'] != 'T/m/s':
        raise ValueError('Please do not modify Siemens slew rate maximum value or units.')
    if sysSiemens['grad_raster_time'] != 10e-6 or sysSiemens['rf_raster_time'] != 1e-6:
        raise ValueError('Please do not modify Siemens gradient/rf raster times.')
    if sysGE['max_grad'] != 31 or sysSiemens['grad_unit'] != 'mT/m':
        raise ValueError('Please do not modify GE gradient maximum value or units.')
    if sysGE['max_slew'] != 130 or sysSiemens['slew_unit'] != 'T/m/s':
        raise ValueError('Please do not modify GE slew rate maximum value or units.')
    if sysGE['grad_raster_time'] != 4e-6 or sysGE['rf_raster_time'] != 4e-6:
        raise ValueError('Please do not modify GE gradient/rf raster times.')

    # Field of View (FOV and resolution)
    acq_params = {}
    acq_params['n'] = 72 # 74 # 72                                        # in-plane  matrix size of reconstructed image
    acq_params['fov'] = 240e-3 # 220e-3 # 240e-3                                  # in-plane fov (m)
    acq_params['nz'] = 54 # 30 # 54                                       # number of reconstructed pixels along z

    if acq_params['n'] <= 0 or acq_params['fov']<= 0 or acq_params['nz'] <= 0:
        raise ValueError('The acquisition parameters can not be zero or negative')

    acq_params['fovz'] = round(acq_params['nz'] * acq_params['fov']*1e2 / acq_params['n']) * 1e-2 # fov along z (m)
    acq_params['dz'] = acq_params['fovz'] / acq_params['nz']    # reconstructed slice thickness (m)
    acq_params['dx'] = acq_params['fov'] / acq_params['n']      # in-plane voxel dimension (m)

    # Slab-selective excitation
    rf_tipdown = {}
    rf_tipdown['flip'] = 10                                     # excitation angle (degrees)
    rf_tipdown['slab_thickness'] = 0.8 * acq_params['fovz']     # slab thickness (m)
    rf_tipdown['tbw'] = 8                                       # time-bandwidth product of SLR pulse
    rf_tipdown['dur'] = 1e-3                                    # RF pulse duration (sec)

    # Fat saturation pulse
    rf_fatsat = {}
    rf_fatsat['flip'] = 50
    # TODO: tbw has to be at least 1.5
    rf_fatsat['tbw'] = 1.5
    rf_fatsat['dur'] = 3e-3

    if rf_tipdown['flip'] < 0 or rf_tipdown['flip'] > 360 or rf_fatsat['flip'] < 0 or rf_fatsat['flip'] > 360:
        raise ValueError('Flip angle has to be between 0-360 degrees')
    if rf_tipdown['dur'] <= 0 or rf_fatsat['dur'] <= 0:
        raise ValueError('RF duration has to be larger than zero')

    # Sequence-specific parameters
    fmri = {}
    fmri['type'] = 'presto'

    if fmri['type'] != 'presto' and fmri['type'] != 'spgr':
        raise ValueError('This type of contrast is not supported. Please select either presto or spgr.')

    fmri['nCyclesSpoil'] = 1 # Spoiler gradient sizes (cycles/voxel). Played on x and z axes.
    # Value of about 1.0-1.5 Gives near-optimal temporal SNR for PRESTO fMRI
    fmri['rf_spoil_seed'] = 150
    fmri['nz_samp'] = acq_params['nz'] #30                                        # kz points sampled per time-frame
    fmri['TR'] = 16.7e-3                                        # sequence TR (sec)
    fmri['dur'] = 5 * 60                                        # total duration of fMRI scan (sec)
    fmri['trVol'] = fmri['nz_samp'] * fmri['TR']                # time to acquire one under-sampled volume (sec)
    fmri['nt'] = 3 # 30                                             # number of time-frames

    # Spiral design parameters
    fmri['nLeafs'] = 3                                          # number of spiral leafs for full k-space sampling
    fmri['dsamp'] = 600                                         # number of samples for the fully sampled core

    if fmri['nt'] < fmri['nLeafs']:
        raise ValueError('The number of frames has to be at least equal to the number of leafs.')

    # Fully- sampled acquisition
    kzFull = np.zeros((acq_params['nz'],))
    for ii in range(1,acq_params['nz']+1):
        kzFull[ii - 1] = ((ii - 1 + 0.5) - acq_params['nz'] / 2) / (acq_params['nz'] / 2)  # scaling is (-1 1)

    if kzFull.min() <= -1 or kzFull.max() >= 1:
        raise ValueError('kzFull has to be scaled between (-1, 1). There is something off.')

    fmri['kzFull'] = kzFull                                     # fully sampled kz sampling pattern
    fmri['nref'] = 0                                            # number of fully sampled frames at the beginning
    fmri['nframes'] = fmri['nref'] * fmri['nLeafs'] + fmri['nt'] # total number of frames

    params = {'sysSiemens': sysSiemens, 'sysGE': sysGE, 'acq_params': acq_params, 'fmri': fmri, 'rf_fatsat': rf_fatsat,
              'rf_tipdown': rf_tipdown}

    return params