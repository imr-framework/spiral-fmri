function seq = getparams
% Get acquisition parameters for 3D B0 scan, and for 3D stack-of-spirals fMRI scan.
%
% In addition, write modules.txt.
%
% (I use the struct name 'seq' here but probably best to call it something else 
% when calling this function to avoid conflict with Pulseq code!)

% local to this function
test = true;         % Create a sequence containing only a few time frames (for testing recon, etc)

%% Global parameters (common to both B0 and fMRI sequence)
seq.fmri.writeKspace = true;   % Write k-space locations for entire scan to a (large) .mat file. Units: cycles/cm

% system limit structs
% NB! If passing 'sys' to writemod.m, 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' can always be a design choice, i.e., it can be at or below the physical system limit.
% Here we will therefore use 'sys' when DESIGNING waveforms, and 'sysGE' when WRITING them to .mod files with writemod.m.
mxs = 13.0;    % max slew [G/cm/msec]. Go easy to minimize PNS.
seq.sys = toppe.systemspecs('maxSlew', mxs, 'slewUnit', 'Gauss/cm/ms', 'maxGrad', 3.1, 'gradUnit', 'Gauss/cm');  % used in design of waveforms
seq.sysGE = toppe.systemspecs('maxSlew', mxs, 'slewUnit', 'Gauss/cm/ms', 'maxGrad', 5, 'gradUnit', 'Gauss/cm');  % needed in writemod() calls only 

% Create modules.txt (common to both B0 and fMRI sequence)
% Entries are tab-separated.
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'readout.mod	0	0	1\n' ...
'tipdown.mod	0	1	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

% FOV and resolution
seq.n = 74;                     % in-plane matrix size of reconstructed image
seq.fov = 22;                   % in-plane fov (cm)
seq.nz = 54;                    % number of reconstructed pixels along z
seq.fovz = 16;                  % fov along z (cm)
seq.dz = seq.fovz/seq.nz;       % reconstructed slice thickness (cm)
seq.dx = seq.fov/seq.n;         % in-plane voxel dimension (cm)

% Slab-selective excitation
seq.rf.flip = 10;                  % excitation angle (degrees)
seq.rf.slabThick = 0.8*seq.fovz;   % cm
seq.rf.tbw = 8;                    % time-bandwidth product of SLR pulse 
seq.rf.dur = 1;                    % RF pulse duration (msec)
seq.rf.ftype = 'ls';               % least-squares SLR pulse (another good option for 3D imaging is 'min')

% fat saturation pulse
seq.fatsat.flip = 50;
%slThick = 1000;     % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse
seq.fatsat.tbw = 1.5;
seq.fatsat.dur = 3;            % pulse duration (msec)

%% Sequence-specific parameters

% Spoiler gradient sizes (cycles/voxel). Played on x and z axes.
% Value of about 1.0-1.5 Gives near-optimal temporal SNR for PRESTO fMRI (see one of my ISMRM abstracts, 2017 I think) 
seq.fmri.nCyclesSpoil = 1;   

seq.fmri.nz_samp = 30;             % Sample this many kz points per time-frame (undersampling factor in kz is 54/30 = 1.8)
seq.fmri.TR = 16.7e-3;             % approximate sequence TR (sec). See toppe.getTRtime()
seq.fmri.dur = 5*60;               % total duration of fMRI scan (sec)
seq.fmri.trVol = seq.fmri.nz_samp*seq.fmri.TR;    % time to acquire one under-sampled image volume (sec)
if test
	seq.fmri.nt = 30;
else
	seq.fmri.nt = 2*round(seq.fmri.dur/seq.fmri.trVol/2);      % number of (undersampled) time-frames
end

% spiral design parameters
seq.fmri.nLeafs = 3;               % Number of spiral rotations (leafs) for full k-space sampling.
seq.fmri.dsamp = 600;              % number of samples for the dense (fully sampled) core. See toppe.utils.spiral.genspivd2()

% fully sampled kz sampling pattern
for ii = 1:seq.nz
	seq.fmri.kzFull(ii) = ((ii-1+0.5)-seq.nz/2)/(seq.nz/2);    % scaling is (-1 1)
end

% undersampled kz sampling pattern
if 0
	% Variable-density (non-Cartesian) kz undersampling. May be useful later.
	Rz = seq.nz/nz_samp;               % kz acceleration factor
	kz = vardenskz(seq.nz,Rz,3.3);    % Fully sampled center with quadratically increasing FOV outside center. Last argument is FOV(kz=0)/FOV(kzmax). 
	a_gz_max = abs((0.5-seq.nz/2)/(seq.nz/2));
	seq.fmri.kzU = kz*a_gz_max;     % scaling is (-1 1)
else
	% Cartesian variable-density kz undersampling
	load zInd;
	seq.fmri.kzU = seq.fmri.kzFull(logical(zInd));
end

seq.fmri.nref = 4;    % number of fully sampled frames acquired at beginning (for, e.g., GRAPPA calibration)

seq.fmri.nframes = seq.fmri.nref*seq.fmri.nLeafs + seq.fmri.nt;   % total number of frames

% For 30 kz platters per frame and rf_spoil_seed=150, ...
% we have mod(nz_samp*rf_spoil_seed,360)=180 which is what we need for improved RF spoiling...
% using our method in Magn Reson Med. 2016 Jun;75(6):2388-93. doi: 10.1002/mrm.25843.
seq.fmri.rf_spoil_seed = 150;  

%% B0 mapping sequence parameters
seq.b0.nCyclesSpoil = 1.5;   % SPGR spoiler gradient size
seq.b0.rf_spoil_seed = 117;  
seq.b0.tdelay = 5;           % (msec) Slow down scan a bit to improve SNR and reduce artifacts from imperfect spoiling.

return;

