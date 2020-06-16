"""
Main script for spiral-fmri reconstruction
Author: Marina Manso Jimeno
Last updated: 06/11/2020
"""
import scipy.io as sio
import numpy as np
import configparser
import os

from tqdm import tqdm

from spiralfmri.Acq.pypulseq.getparams import getparams
from spiralfmri.Recon.Python.utils.recon_opts import data_check_file_inputs, data_check_getparams_input, rfphs_compensate, spiral_alignment
from spiralfmri.Recon.Python.recon_SoS import recon_SoS
from spiralfmri.Recon.Python.utils.display_timepoint import disp_timepoint

# Read settings.ini configuration file
path_settings = '../../settings.ini'
config = configparser.ConfigParser()
config.read(path_settings)
config_data = config['RECON']

##
# Input raw data from scanner + k-space trajectory coordinates
##
dat_file = sio.loadmat(config_data['path_rawdat_file'])['ksp_dat'] # shape Npoints x kz x Nframes x Nchannels
# TODO: Modify saving of ktraj during acq and dimensions
ktraj = sio.loadmat(config_data['path_ktraj_file'])['ktraj']
ktraj = np.moveaxis(ktraj, -1, 0)
ktraj = np.moveaxis(ktraj, 1, 2)
data_check_file_inputs(raw_dat=dat_file, ktraj=ktraj)

# Acquisition parameters
##
seq_params = getparams()
data_check_getparams_input(raw_dat=dat_file, params=seq_params)

sys_acq = seq_params['sysSiemens']
sys_gen = seq_params['sysGE']
acq_params = seq_params['acq_params']
fmri = seq_params['fmri']
rf_fatsat = seq_params['rf_fatsat']
rf_tipdown = seq_params['rf_tipdown']

kz = acq_params['nz']
nframes = fmri['nframes']
rfphs = fmri['rf_spoil_seed']
nleafs = fmri['nLeafs']
acqWin = range(85,929)

##
# RF phase compensation
##
dat_file_rfphs = rfphs_compensate(dat_file, rfphs)

##
# Align the spiral shots (sliding window reconstruction)
##
dat_file_2recon = spiral_alignment(dat_file_rfphs, nleafs)
dat_file_2recon = dat_file_2recon[acqWin, :, :, :, :]

##
# Arrange the k-space trajectory
##
ktraj = ktraj[acqWin, 0, 0:nleafs, 0:2] # cycles/m
ktraj = np.squeeze(np.apply_along_axis(lambda args: [complex(*args)], 2, ktraj))


##
# Recon
##
N = acq_params['n']
ntimepoints = dat_file_2recon.shape[3]
recon_vol = np.zeros((N, N, kz, ntimepoints))
for tp in tqdm(range(ntimepoints)):
    recon_vol[:,:,:,tp] = recon_SoS(dat_file_2recon[:,:,:,tp,:], ktraj, N)
    disp_timepoint(recon_vol[:,:,:,tp])

# Save the result
np.save(os.path.join(config_data['path_save_recon'], 'recon_vol.npy'), recon_vol)

