# Copyright of the Board of Trustees of Columbia University in the City of New York

"""
PRESTO stack of spirals pypulseq version
Author: Marina Manso Jimeno
Last Updated: 05/07/2020
"""
# Import Pypulseq methods
from pypulseq.opts import Opts
from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import make_adc
from pypulseq.make_arbitrary_rf import make_arbitrary_rf
from pypulseq.make_arbitrary_grad import make_arbitrary_grad

# Import own methods
from getparams import getparams
from g2k import g2k
#from genspivd2 import genspiralvd
from interpolation import grad_interpolate
from interpolation import rf_interpolate
from archimedian import archimedian_spiral
from make_slr_rf import make_slr_rf
from makeSystemlength import makeSystemlength
from make_balanced import make_balanced
from make_crusher import make_crusher
from trapwave2 import trapwave2

# Other imports
from math import sqrt, pi
from cmath import exp as cexp
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


###
# 1. Load the acquisition parameters
###
seq_params = getparams()
sys_acq = seq_params['sysSiemens']
sys_gen = seq_params['sysGE']
acq_params = seq_params['acq_params']
fmri = seq_params['fmri']
rf_fatsat = seq_params['rf_fatsat']
rf_tipdown = seq_params['rf_tipdown']

sys = sys_acq
system = Opts(max_grad=sys['max_grad'], grad_unit=sys['grad_unit'], max_slew=sys['max_slew'], slew_unit=sys['slew_unit'], grad_raster_time=sys['grad_raster_time'], rf_raster_time=sys['rf_raster_time'])
seq = Sequence(system)

###
# 2. Generate the RF and gradient waveforms
###
'''
Fat Saturation Module
'''
slThick = 10 # dummy value [m]
rf_fs_wav, _ = make_slr_rf(rf_fatsat['flip'], slThick, rf_fatsat['tbw'], rf_fatsat['dur'], 0, sys_gen)
rf_fs_wav = rf_interpolate(rf_fs_wav, sys_acq['rf_raster_time'])
'''
Slab selective excitation module and PRESTO gradients (Tip down module)
'''
# RF pulse
nCyclesSpoil = 0
rf_td_wav, g_td_wav = make_slr_rf(rf_tipdown['flip'], rf_tipdown['slab_thickness'], rf_tipdown['tbw'], rf_tipdown['dur'], nCyclesSpoil, sys_gen)

# Spoiler gradients
gspoil1 = make_crusher(fmri['nCyclesSpoil'], acq_params['dz'], 0, sys_gen['max_grad'], 0.9 * sys_gen['max_slew'] / sqrt(2), sys_gen)
gspoil2 = make_crusher(2*fmri['nCyclesSpoil'], acq_params['dz'], 0, sys_gen['max_grad'], 0.7 * sys_gen['max_slew'] / sqrt(2), sys_gen)
# TODO: Interpolate grad
if fmri['type'] == 'spgr':
    gspoil2 = gspoil1
    gspoil1 = []

# tipdown module
rf_td_wav = np.vstack((0 * gspoil2, np.zeros((2,1)), rf_td_wav, 0 * gspoil1))
gx_td_wav = np.vstack((gspoil2, np.zeros((2,1)), 0 * g_td_wav, -gspoil1))
gz_td_wav = np.vstack((gspoil2, np.zeros((2,1)), g_td_wav, -gspoil1))

rf_td_wav = makeSystemlength(rf_td_wav, sys_gen['rf_raster_time'])
gx_td_wav = makeSystemlength(gx_td_wav, sys_gen['grad_raster_time'])
gz_td_wav = makeSystemlength(gz_td_wav, sys_gen['grad_raster_time'])

rf_td_wav = rf_interpolate(rf_td_wav, sys_acq['rf_raster_time'])
gx_td_wav = grad_interpolate(gx_td_wav)
gz_td_wav = grad_interpolate(gz_td_wav)

'''
Balanced stack-of-spirals (Readout module)
'''
# Two concatenated archemidian spirals with two velocities
# g_ro_wav, ksp, _, _ = genspiralvd(acq_params['fov'], acq_params['n'], 3,  0.99 * system.max_grad, 0.99 * system.max_slew, sys['grad_raster_time'], fmri['dsamp'])

# Archimedian spiral (ISMRM version) in Python
rmax = acq_params['n'] / (2 * acq_params['fov']) # max k-space radius
k, g_ro_wav = archimedian_spiral(0.99 * sys_gen['max_slew'], 0.99 * sys_gen['max_grad'], sys_gen['grad_raster_time'], fmri['nLeafs'], acq_params['fov'], rmax)

# add a couple of zeros to make sure k=0 is sampled
g_ro_wav = np.vstack((np.zeros((2, 1)), g_ro_wav.reshape(g_ro_wav.shape[0], 1)))  # Hz/m

# make balanced and same length
gx_ro_wav = make_balanced(g_ro_wav.real, sys_acq['max_grad'], sys_acq['max_slew']/sqrt(2), sys_gen)
gy_ro_wav = make_balanced(g_ro_wav.imag, sys_acq['max_grad'], sys_acq['max_slew']/sqrt(2), sys_gen)
n = max(len(gx_ro_wav), len(gy_ro_wav))
gx_ro_wav = np.concatenate((gx_ro_wav, np.zeros((n - len(gx_ro_wav),1))))
gy_ro_wav = np.concatenate((gy_ro_wav, np.zeros((n - len(gy_ro_wav),1))))

# make it spiral in
gx_ro_wav = np.flipud(gx_ro_wav)
gy_ro_wav = np.flipud(gy_ro_wav)

# add partition (kz) enconding trapezoids
gzamp = (1/sys_gen['grad_raster_time']) / (system.gamma * acq_params['fovz']) * 1e3 # mT/m
zarea = gzamp * acq_params['nz'] * sys_gen['grad_raster_time'] # mT/m*s
gpe = - trapwave2(zarea / 2, sys_acq['max_grad'], sys_acq['max_slew'], sys_gen['grad_raster_time']) * 1e-3 * system.gamma # Hz/m
gx_ro_wav_orig = makeSystemlength(np.vstack((0 * gpe.reshape(-1,1), np.zeros((2,1)), gx_ro_wav, 0*gpe.reshape(-1,1))), sys_gen['grad_raster_time'])#gx_ro_wav.T
gy_ro_wav_orig = makeSystemlength(np.vstack((0 * gpe.reshape(-1,1), np.zeros((2,1)), gy_ro_wav, 0*gpe.reshape(-1,1))), sys_gen['grad_raster_time'])#gy_ro_wav.T
gz_ro_wav_orig = makeSystemlength(np.vstack((gpe.reshape(-1,1), np.zeros((2,1)), 0 * gx_ro_wav, -gpe.reshape(-1,1))), sys_gen['grad_raster_time'])

gx_ro_wav_orig = grad_interpolate(gx_ro_wav_orig)
gy_ro_wav_orig = grad_interpolate(gy_ro_wav_orig)
gz_ro_wav_orig = grad_interpolate(gz_ro_wav_orig)

# ADC
adc = make_adc(num_samples=max(gx_ro_wav_orig.shape), dwell=system.grad_raster_time)

###
# 3. Sequence loop
###
rfphs = 0  # rad
rf_spoil_seed_cnt = 0

ktraj_full = []
for iframe in range(1, fmri['nframes'] + 1):
    print(str(iframe) + ' of ' + str(fmri['nframes']))

    # set kz sampling pattern for this frame
    kz1 = fmri['kzFull']  # numel(kz) = nz

    ktraj_frame = []
    for iz in range(1, len(kz1)+1):

        # fat sat block
        rf_fs, _ = make_arbitrary_rf(flip_angle= rf_fatsat['flip'] * 180 / pi, system= system, signal=np.squeeze(rf_fs_wav), freq_offset= -440, phase_offset= rfphs)#, time_bw_product=1.5)

        seq.add_block(rf_fs)

        # Slab excitation + PRESTO gradients block
        rf_td, _ = make_arbitrary_rf(flip_angle= rf_tipdown['flip'] * 180 / pi, system= system, signal= np.squeeze(rf_td_wav), phase_offset= rfphs)#,time_bw_product=8)

        gx_td = make_arbitrary_grad(channel= 'x', system= system, waveform= np.squeeze(gx_td_wav))

        gz_td = make_arbitrary_grad(channel= 'z', system= system, waveform= np.squeeze(gz_td_wav))

        seq.add_block(rf_td, gx_td, gz_td)

        # Read out block
        slice = iz
        ileaf = (((iframe % fmri['nLeafs']) + iz) % fmri['nLeafs']) + 1  # rotate leaf every frame and every kz platter
        phi = 2 * pi * (ileaf - 1) / fmri['nLeafs'] # leaf rotation angle in radians

        # Rotate the spiral arm if necessary
        if phi != 0:
            g = np.zeros(gx_ro_wav_orig.shape, dtype=complex)
            for i in range(gx_ro_wav_orig.shape[0]):
                g[i] = complex(gx_ro_wav_orig[i], gy_ro_wav_orig[i])
                g[i] = g[i] * cexp(phi * 1j)
            gx_ro_wav = np.real(g)
            gy_ro_wav = np.imag(g)

        else:
            gx_ro_wav = gx_ro_wav_orig
            gy_ro_wav = gy_ro_wav_orig

        gx_ro = make_arbitrary_grad(channel= 'x', system= system, waveform= np.squeeze(gx_ro_wav))

        gy_ro = make_arbitrary_grad(channel= 'y', system= system, waveform= np.squeeze(gy_ro_wav))

        gz_ro_wav = gz_ro_wav_orig * kz1[slice-1]
        gz_ro = make_arbitrary_grad(channel= 'z', system= system, waveform= np.squeeze(gz_ro_wav))

        seq.add_block(gx_ro, gy_ro, gz_ro, adc)

        # get the k-space trajectory
        kx, ky, ks_z = g2k(gx_ro.waveform, gy_ro.waveform, gz_ro.waveform)
        ktraj_frame.append(np.array([kx, ky, ks_z]))

        # update rf phase (RF spoiling)
        rfphs = rfphs + (fmri['rf_spoil_seed'] / 180 * pi) * rf_spoil_seed_cnt  # rad
        rf_spoil_seed_cnt += 1

    ktraj_full.append(ktraj_frame)

###
# 4. Plot the resulting sequence
###
seq.plot(time_range=(0,0.019))

###
# 5. Write the .seq file
###
seq.write('spiral_' + fmri['type'] + '.seq')
sio.savemat('ktraj.mat', {'ktraj_full': ktraj_full})



