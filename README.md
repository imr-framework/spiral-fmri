# spiral-fMRI
Pulseq implementation of the existing TOPPE spriral PRESTO sequence for fMRI written in MATLAB. Visit [their GitLab repository](https://gitlab.com/fMRI/toppe-sequences/stack-of-spirals-presto-bold-fmri.git). 

## Why?
The availability of open-source software for rapid prototyping of MRI pulse sequences has increased in the last few years. Two examples are the [TOPPE](https://www.ncbi.nlm.nih.gov/pubmed/29096052) and [Pulseq](https://www.ncbi.nlm.nih.gov/pubmed/27271292) frameworks.

Both platforms facilitate considerably the implementation of non-standard sequences in commercial scanners such as GE, Siemens and Bruker systems. However,the same sequence can have different outcomes when using different vendors due to the subtle differences in system specifications and sequence parameters. 

This project aims to test Pulseq's performance in scanners from two different vendors, specifically GE and Siemens. For that purpose a stack-of-spirals [PRESTO sequence for fMRI](https://www.ncbi.nlm.nih.gov/pubmed/22245350) has been "translated" from TOPPE (MATLAB written) to [pyPulseq](https://github.com/imr-framework/pypulseq) and made available for anyone who wants to give it a try. 

Enjoy!

Please visit the hyperlinks if you want to learn more.

## How?

**Please note:** You need to install Pulseq interpreter on your system before running any Pulseq sequence.

In **PRESTO_pulseq**:
1. Change the system limits if you need to. 
    ```python
    gamma = 42576000  # in Hz/T  %Determined from Pulseq - do not change

    kwargs_for_opts = {'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s', 'grad_dead_time': 10e-6}
    system = Opts(kwargs_for_opts)
    seq = Sequence(system)
    ```
    
    For example, gradient and slew rate limits are not the same for GE and Siemens.
    
    ![GEvsSiemens system limits](/images/rf_g_limits.png)
   
   
2. Run the code . It will give you the .seq file that will be played on the scanner.
    ```Python
    seq.write("spiral_PRESTO.seq")
    ```
3. Go to the MR room and play it!

## More about the code

### 1. RF and gradient Interpolation: 
As noticed before, RF and gradient raster times are not the same across vendors. Therefore, waveforms have to be interpolated.
  ```Python
  # Import waveforms
  gx_ro_wav = sio.loadmat('gx1_readout.mat')['gx1'] * gamma / 100  # from G/cm to Hz/m
  gy_ro_wav = sio.loadmat('gy1_readout.mat')['gy1'] * gamma / 100
  gz_ro_wav = sio.loadmat('gz1_readout.mat')['gz1'] * gamma / 100

  # Interpolate and reshape
  gx_ro_wav_orig = grad_interpolate(gx_ro_wav, seq.grad_raster_time).reshape((1, -1))
  gy_ro_wav_orig = grad_interpolate(gy_ro_wav, seq.grad_raster_time).reshape((1, -1))
  gz_ro_wav_orig = grad_interpolate(gz_ro_wav, seq.grad_raster_time).reshape((1, -1
  ```
  
  ![RF interpolation](/images/rf_interp.png)
  
  ![Gradient interpolation](/images/g_interp.png)
  
  ### 2. Spiral Rotation
  Spiral waveforms have to be rotated in order to get an uniform k-space coverage. 
  
```Python
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
```
      
  ![3 Superimposed spiral arms](/images/spiral_rot.png)
  ### 3. k-space trajectory
  Getting the k-space trajectory is important to get a quick sense of what you are running on the scanner and for later reconstruction.
  
```Python
# get the k-space trajectory
kx, ky, ks_z = g2k(gx_ro.waveform, gy_ro.waveform, gz_ro.waveform)

ktraj_frame.append(np.array([kx, ky, ks_z]))
```
   
   ![k-space trajectory](/images/1_9_19_29_39_49.png)
        

