# Stack-of-spirals PRESTO fMRI for GE (TOPPE)

## PREREQUISITES

### TOPPE Matlab toolbox
See https://toppemri.github.io/
```
$ git clone git@github.com:toppeMRI/toppe.git
```

Add the 'toppe' folder to your Matlab path.
```matlab
>> addpath('./toppe/');
```
### 3D GRAPPA recon Matlab toolbox

TODO


## LOW-RESOLUTION SPOILED 3D GRE (SPGR/FLASH) 

Use to obtain:
* B0 map
* Coil sensitivity maps (e.g., using BART, see ./recon/sense/sensmaps-espirit/main.m) TODO
	
To create sequence files (Matlab):
```matlab
>> params = getAcParams('b0map'); 
>> createGEfiles(params);
```

This creates ./tar/scan,b0.tar. Untar this file in /usr/g/bin/ on scanner and scan with toppev3.e.


#### To reconstruct

TODO


## STACK-OF-SPIRALS 3D PRESTO fMRI

Dynamic T2\*-weighted PRESTO sequence
* 3.33 mm isotropic resolution
* Fully-sampled TR is about 18e-3x(3x60+nskip) = ~3sec
* Undersampling factor is R=5.4 (3 in-plane x 1.8 in kz) => TR is about 1/2 sec

To create sequence files (Matlab):
```matlab
>> params = getAcParams('fmri'); 
>> createGEfiles(params);
```

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
>> params = getAcParams('fmri'); 
>> ims = recon(d,ksp,params.fov,params.imSize);
```


