# Stack-of-spirals PRESTO fMRI for GE (TOPPE)

## PREREQUISITES

### TOPPE Matlab toolbox
See https://toppemri.github.io/
```
$ git clone git@github.com:toppeMRI/toppe.git
```

It's a good idea to update your local copy periodically:
```
$ git pull
```

Add the 'toppe' folder to your Matlab path.
```matlab
>> addpath('./toppe/');
```
### 3D GRAPPA Matlab toolbox

TODO


## MAIN SCRIPT

```matlab
>> params = getparams;
>> main(params);
```

This creates two sequences:
* ./tar/scan,b0.tar
  * Use to obtain B0 map and coil sensitivity maps
* ./tar/scan,fmri.tar

Untar each file (in turn) in /usr/g/bin/ on scanner host and scan with toppev3.e.

Notes:
* Modules are designed so that **predicted PNS is < 80% of stimulation threshold** for all modules.
To view predicted PNS, see toppe.pns, or do
```matlab
>> toppe.plotmod('all');
```

## LOW-RESOLUTION SPOILED 3D GRE (SPGR/FLASH) 

TODO


## STACK-OF-SPIRALS 3D PRESTO fMRI

Dynamic T2\*-weighted PRESTO sequence
* 3.33 mm isotropic resolution
* Fully-sampled TR is about 18e-3x(3x60+nskip) = ~3sec
* Undersampling factor is R=5.4 (3 in-plane x 1.8 in kz) => TR is about 1/2 sec

#### Data pre-processing
Apply temporal filter to raw k-space data BEFORE reconstructing, to remove ghosting, see Magn Reson Med. 2016 Jun;75(6):2388-93. doi: 10.1002/mrm.25843.
This trick relies on using the same k-space trajectory for each (undersampled) temporal frame.
```matlab
>> d = toppe.utils.loadpfile(pfile);
>> d = tfilt(d);
```

#### Image reconstruction
```matlab
>> ksp = getspiralkspace;
>> params = getparams;
>> imSize = [params.n params.n params.nz];
>> ims = recon(d, ksp, params.fov, imSize);
```


