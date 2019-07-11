"""
PRESTO stack of spirals pulseq version
Author: Marina Manso Jimeno
Last Updated: 06/12/2019
"""
"""
Import general methods, pypulseq methods and functions
"""
from itertools import compress
from math import pi
from cmath import exp

import numpy as np
import scipy.io as sio

from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import makeadc
from pypulseq.makearbitrary_grad import makearbitrary_grad
from pypulseq.makearbitrary_rf import make_arbitrary_rf
from pypulseq.opts import Opts

from interpolation import rf_interpolate,grad_interpolate
from g2k import g2k


"""
Testing mode
"""
test = True  # Create a sequence containing only a few time frames (for testing recon, etc)

"""
Gradient limits definition for maximum slew rate and maximum amplitude
"""

gamma = 42576000  # in Hz/T  %Determined from Pulseq - do not change

kwargs_for_opts = {'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s', 'grad_dead_time': 10e-6}
system = Opts(kwargs_for_opts)
seq = Sequence(system)

"""
Acquisition parameters definition
Currently RF and gradients waveforms are imported from the TOPPE Matlab code version

fov = 240e-3 # in-plane fov
N = 72 # in-plane matrix size of reconstructed image
Nz = 54 # number of reconstructed pixels along z. (Excitation is 14 cm, readout FOV in z is 54*3.33mm = 180mm) #TODO: Nz
fovz = 180e-3 # fov along z
slice_thickness = fovz / Nz # reconstructed slice thickness

n_shots = 3 # number of spiral rotations (leafs) for full k-space sampling
alpha = 8 # oversampling factor
"""

Nz = 54
n_shots = 3
rf_td_flip = 10 * 180 / pi  # rad
rf_fs_flip = 50 * 180 / pi  # rad

"""
Import the different RF and gradient waveforms
RF waveforms are converted from [G] to [Hz] and interpolated for 1 microsecond raster time. 
Gradient waveform are converted from [G/cm] to [Hz/m] and interpolated for 10 microseconds raster time
All waveforms are also reshaped as (1,x) to fulfill the seq.add_block method requirements
"""

"""
Read-out block: balanced stack-of-spirals  
"""
# Import waveforms
gx_ro_wav = sio.loadmat('gx1_readout.mat')['gx1'] * gamma / 100  # from G/cm to Hz/m
gy_ro_wav = sio.loadmat('gy1_readout.mat')['gy1'] * gamma / 100
gz_ro_wav = sio.loadmat('gz1_readout.mat')['gz1'] * gamma / 100

# Interpolate and reshape
gx_ro_wav_orig = grad_interpolate(gx_ro_wav, seq.grad_raster_time).reshape((1, -1))
gy_ro_wav_orig = grad_interpolate(gy_ro_wav, seq.grad_raster_time).reshape((1, -1))
gz_ro_wav_orig = grad_interpolate(gz_ro_wav, seq.grad_raster_time).reshape((1, -1))

# ADC
kwargs_for_adc = {"num_samples": max(gx_ro_wav_orig.shape), "dwell": system.grad_raster_time}
adc = makeadc(kwargs_for_adc)

"""
Slab-selective excitation (tip-down) + PRESTO spoiler gradients
"""
# Import waveforms
rf_td_wav = sio.loadmat('rf_tipdown.mat')['rf'] * gamma / 10000  # from G to Hz
gx_td_wav = sio.loadmat('gx_tipdown.mat')['gx'] * gamma / 100  # from G/cm to Hz/m
gz_td_wav = sio.loadmat('gz_tipdown.mat')['gz'] * gamma / 100

#Interpolate and reshape
rf_td_wav = rf_interpolate(rf_td_wav, seq.rf_raster_time).reshape((1, -1))
gx_td_wav = grad_interpolate(gx_td_wav, seq.grad_raster_time).reshape((1, -1))
gz_td_wav = grad_interpolate(gz_td_wav, seq.grad_raster_time).reshape((1, -1))

"""
Fat saturation pulse
"""
# Import waveforms
rf_fs_wav = sio.loadmat('rf_fatsat.mat')['rf'] * gamma / 10000  # from G to Hz

#Interpolate and reshape
rf_fs_wav = rf_interpolate(rf_fs_wav, seq.rf_raster_time).reshape((1, -1))

#Frequency offset to null fat signal
fs_freq_offset = -440 # Hz

"""
Sequence creation
"""
nz_samp = 30  # kz points sampled per time-frame (undersampling factor 54/30 = 1.8)
dur = 5 * 60  # total duration of the fMRI scan (sec)

# Undersampling factors

Rxy = 3
Rz = Nz / nz_samp

trVol = 2.9  # time for fully sampled image volume (sec) (approximate)

if test:
    nt = 3 * 10
else:
    nt = 2 * round(Rxy * Rz * dur / trVol / 2)  # number of undersampled time-frames

# fully sampled kz sampling pattern
kzFull = []
for i in range(1, Nz + 1):
    kzFull.append(((i - 1 + 0.5) - Nz / 2) / (Nz / 2))  # scaling is (-1,1)

# undersampled kz sampling pattern
if 0:
    pass
    # variable density (non-cartesian) kz undersampling. May be useful later
else:
    # cartesian variable-density kz undersampling
    zInd = sio.loadmat('zInd.mat')['zInd']
    zInd = (zInd == 1).tolist()[0]
    kzU = list(compress(kzFull, zInd))

ndisdaq = 0 # was 10 before
nref = 4 * n_shots  # 4 is the number of fully sampled frames acquired at beginning

rfphs = 0  # rad
rfphsLast = rfphs
daqphs = 0
rf_spoil_seed_cnt = 0
rf_spoil_seed = 150  # % For 30 kz platters per frame and rf_spoil_seed=150, we have mod(nz_samp*rf_spoil_seed,360)=180 which is what we need for improved RF spoiling using our method in Magn Reson Med. 2016 Jun;75(6):2388-93. doi: 10.1002/mrm.25843.

# loop over undersampled time frames
dabon = 1
daboff = 0
ktraj_full = []
for iframe in range(-ndisdaq + 1, nref + 1):#range(-ndisdaq + 1, nref + nt + 1):
    print(str(iframe + ndisdaq) + ' of ' + str(nt + ndisdaq + nref))

    # set kz sampling pattern for this frame
    if iframe < (nref + 1):
        kz = kzFull  # numel(kz) = nz
    else:
        kz = kzU  # numel(kz) = nz_samp

    # set  'view' data storage index
    if iframe < 1:
        dabmode = 'off'
        view = 1
    else:
        dabmode = 'on'
        view = iframe

    ktraj_frame =[]
    for iz in range(1, len(kz) + 1):

        # fat sat block

        kwargs_for_arbRF = {'flip_angle': rf_fs_flip, 'system': system, 'signal': rf_fs_wav, 'freq_offset': fs_freq_offset, 'phase_offset': rfphs}
        rf_fs, _ = make_arbitrary_rf(kwargs_for_arbRF)

        seq.add_block(rf_fs)

        # Slab excitation + PRESTO gradients block
        kwargs_for_arbRF = {'flip_angle': rf_td_flip, 'system': system, 'signal': rf_td_wav, 'phase_offset': rfphs}
        rf_td, _ = make_arbitrary_rf(kwargs_for_arbRF)

        kwargs_for_arbG = {'channel': 'x', 'system': system, 'waveform': gx_td_wav}
        gx_td = makearbitrary_grad(kwargs_for_arbG)

        kwargs_for_arbG = {'channel': 'z', 'system': system, 'waveform': gz_td_wav}
        gz_td = makearbitrary_grad(kwargs_for_arbG)

        seq.add_block(rf_td, gx_td, gz_td)

        # Read out block
        slice = iz
        if iframe < (nref + 1):
            ileaf = (((iframe % n_shots) + iz) % n_shots) + 1  # rotate leaf every frame and every kz platter
        else:
            ileaf = (iz % n_shots) + 1  # rotate every kz platter

        phi = 2 * pi * (ileaf - 1) / n_shots # leaf rotation angle in radians

        # Rotate the spiral arm if necessary
        if phi != 0:
            g = []
            for i in range(len(gx_ro_wav_orig[0])):
                g.append(complex(gx_ro_wav_orig[0, i], gy_ro_wav_orig[0, i]))
                g[i] = g[i] * exp(phi * 1j)
            gx_ro_wav = np.real(g)
            gy_ro_wav = np.imag(g)

        else:
            gx_ro_wav = gx_ro_wav_orig
            gy_ro_wav = gy_ro_wav_orig


        kwargs_for_arbG = {'channel': 'x', 'system': system, 'waveform': gx_ro_wav}
        gx_ro = makearbitrary_grad(kwargs_for_arbG)

        kwargs_for_arbG = {'channel': 'y', 'system': system, 'waveform': gy_ro_wav}
        gy_ro = makearbitrary_grad(kwargs_for_arbG)

        gz_ro_wav = gz_ro_wav_orig * kz[slice-1]
        kwargs_for_arbG = {'channel': 'z', 'system': system, 'waveform': gz_ro_wav}
        gz_ro = makearbitrary_grad(kwargs_for_arbG)

        seq.add_block(gx_ro, gy_ro, gz_ro, adc)

        # get the k-space trajectory
        kx, ky, ks_z = g2k(gx_ro.waveform, gy_ro.waveform, gz_ro.waveform)

        ktraj_frame.append(np.array([kx, ky, ks_z]))

        # update rf phase (RF spoiling)

        rphsLast = rfphs
        rfphs = rfphs + (rf_spoil_seed / 180 * pi) * rf_spoil_seed_cnt  # rad
        rf_spoil_seed_cnt += 1


    ktraj_full.append(ktraj_frame)
"""
Plot the resulting sequence
"""
seq.plot( time_range=(0, 0.019))

"""
Write the .seq file
"""
seq.write("spiral_PRESTO.seq")