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
seq.writeKspace = false;   % Write k-space locations for entire scan to a (large) .mat file. Units: cycles/cm

% system limit structs
% NB! If passing 'sys' to writemod.m, 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' can always be a design choice, i.e., it can be at or below the physical system limit.
% Here we will therefore use 'sys' when DESIGNING waveforms, and 'sysGE' when WRITING them to .mod files with writemod.m.
mxs = 10.0;    % max slew [G/cm/msec]. Go easy to minimize PNS.
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
seq.n = 72;                     % in-plane matrix size of reconstructed image
seq.fov = 24;                   % in-plane fov (cm)
seq.nz = 54;                    % number of reconstructed pixels along z
seq.fovz = 18;                  % fov along z (cm)
seq.dz = seq.fovz/seq.nz;       % reconstructed slice thickness (cm)
seq.dx = seq.fov/seq.n;         % in-plane voxel dimension (cm)

% spiral readout
seq.nleafs = 3;                 % Number of spiral rotations (leafs) for full k-space sampling.

% Slab-selective excitation
seq.rf.flip = 10;                  % excitation angle (degrees)
seq.rf.slabThick = 0.8*seq.fovz;   % cm
seq.rf.tbw = 8;                    % time-bandwidth product of SLR pulse 
seq.rf.dur = 1;                    % RF pulse duration (msec)
seq.rf.nCyclesSpoil = 0;           % make it balanced initially -- gradient pre/rephasers will then be modified
seq.rf.ftype = 'ls';               % least-squares SLR pulse (another good option for 3D imaging is 'min')

% Spoiler gradient sizes (cycles/voxel). Played on x and z axes.
seq.fmri.nCyclesSpoil = 1;   % Gives near-optimal temporal SNR for PRESTO fMRI, and not so large (see one of my ISMRM abstracts, 2017 I think) 
seq.b0.nCyclesSpoil = 1.5;   % SPGR spoiler gradient

% fat saturation pulse
seq.fatsat.flip = 50;
%slThick = 1000;     % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse
seq.fatsat.tbw = 1.5;
seq.fatsat.dur = 3;            % pulse duration (msec)

%% fMRI sequence parameters

seq.fmri.nz_samp = 30;             % Sample this many kz points per time-frame (undersampling factor in kz is 54/30 = 1.8)
seq.fmri.TR = 16.7e-3;             % approximate sequence TR (sec). See toppe.getTRtime()
seq.fmri.dur = 5*60;               % total duration of fMRI scan (sec)
seq.fmri.trVol = nz_samp*seqTR;    % time to acquire one under-sampled image volume (sec)
if test
	seq.nt = 30;
else
	seq.nt = 2*round(seq.fmri.dur/seq.fmri.trVol/2);      % number of (undersampled) time-frames
end

% fully sampled kz sampling pattern
for ii = 1:seq.nz
	seq.kzFull(ii) = ((ii-1+0.5)-seq.nz/2)/(seq.nz/2);    % scaling is (-1 1)
end

% undersampled kz sampling pattern
if 0
	% Variable-density (non-Cartesian) kz undersampling. May be useful later.
	Rz = seq.nz/nz_samp;               % kz acceleration factor
	kz = vardenskz(seq.nz,Rz,3.3);    % Fully sampled center with quadratically increasing FOV outside center. Last argument is FOV(kz=0)/FOV(kzmax). 
	a_gz_max = abs((0.5-seq.nz/2)/(seq.nz/2));
	kzU = kz*a_gz_max;     % scaling is (-1 1)
else
	% Cartesian variable-density kz undersampling
	load zInd;
	seq.kzU = kzFull(logical(zInd));
end

seq.fmri.ndisdaq = 2;        % number of (fully sampled) dummy TRs to reach steady state
seq.fmri.nref    = 4;        % number of fully sampled frames acquired at beginning (for, e.g., GRAPPA calibration)

seq.fmri.nframes = seq.fmri.nref*seq.nleafs + seq.nt;

% For 30 kz platters per frame and rf_spoil_seed=150, ...
% we have mod(nz_samp*rf_spoil_seed,360)=180 which is what we need for improved RF spoiling...
% using our method in Magn Reson Med. 2016 Jun;75(6):2388-93. doi: 10.1002/mrm.25843.
seq.fmri.rf_spoil_seed = 150;  

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write 3D B0 mapping sequence                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create readout module (readout.mod). Matrix size must match the fMRI acquisition.
% Remember: makegre creates y and z phase-encodes based on isotropic resolution
% Again, use 'sys' in design, but 'sysGE' in writemod().
[gx,gy,gz] = toppe.utils.makegre(seq.fov, seq.n, seq.dx, 'system', seq.sys, 'ncycles', 1);
[~,~,~,~,~,hdrints,hdrfloats] = toppe.readmod('readout.mod');
npre = hdrints(1);       % number of samples before readout plateau
npixro = hdrints(2);     % number of samples during readout plateau
oprbw = hdrfloats(20);   % readout bandwidth (max +/ 125 kHz)
toppe.writemod('gx', gx, 'gy', gy, 'gz', gz, 'ofname', 'readout.mod', 'hdrints', [npre npixro], 'hdrfloats', [oprbw], ...
               'desc', '3D spin-warp readout', 'system', seq.sysGE);  % NB! Pass 'sysGE' here, not 'sys'!

%% Create slab-selective excitation module (tipdown.mod)
seqb0.rf.flip = 4;                            % excitation angle (degrees)
seqb0.rf.slabThick = seq.rf.slabThick*1.1;    % cm. Excite a bit more than the fMRI slab to ensure that edge slices don't drop out.
seqb0.rf.tbw = 8;                             % time-bandwidth product of SLR pulse 
seqb0.rf.dur = 1;                             % RF pulse duration (msec)
seqb0.rf.ncyclesspoil = 1e-3;                 % small (but non-zero) so makeslr doesn't put a spoiler (or balancer) in front
seqb0.rf.ftype = 'ls';
[rf,gex] = toppe.utils.rf.makeslr(seqb0.rf.flip, seqb0.rf.slabThick, seqb0.rf.tbw, seqb0.rf.dur, seqb0.rf.ncyclesspoil, ...
             'ftype', seqb0.rf.ftype, 'system', seq.sys, 'ofname', 'tipdown.mod');

seqb0.nCyclesSpoil = 2;
gspoil = toppe.utils.makecrusher(seqb0.nCyclesSpoil, seq.dz, 0, seq.sys.maxSlew, seq.sys.maxGrad);

rf = toppe.utils.makeGElength([0*gspoil(:); zeros(2,1); rf(:)]);
gx = toppe.utils.makeGElength([1*gspoil(:); zeros(2,1); 0*gex(:)]);
gz = toppe.utils.makeGElength([1*gspoil(:); zeros(2,1); gex(:)]);

toppe.writemod('rf', rf, 'gx', gx, 'gz', gz, 'ofname', 'tipdown.mod', ...
               'desc', 'RF slab excitation', 'system', seq.sysGE);  % NB! Pass 'sysGE' here, not 'sys'!


%% Create scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;
textra = 140;          % (msec) slow down sequence to reduce SAR

nx = seq.n;
ny = seq.n;
nz = seq.nz;

tdelay = 3;           % (msec) Slow down scan a bit to improve SNR and reduce artifacts from imperfect spoiling.

fprintf('Writing scanloop.txt for B0 sequence\n');
toppe.write2loop('setup');
for iz = -0:nz           % iz < 1 are discarded acquisitions (to reach steady state)
	fprintf([repmat('\b',1,20) sprintf('%d of %d', iz, nz)]);

	for iy = 1:ny
		for iim = 1:2
			switch iim
				case 1
					textra1 = 0;
					textra2 = 2.3 + tdelay;     % msec
				case 2
					textra1 = 2.3;     % msec
					textra2 = 0 + tdelay;
			end

         if iz > 0
            a_gy = ((iy-1+0.5)-ny/2)/(ny/2);    % y phase-encode amplitude, scaled to (-1,1) range
        		a_gz = ((iz-1+0.5)-nz/2)/(nz/2);    % z phase-encode amplitude, scaled to (-1,1) range
				dabmode = 'on';
         else
            a_gy = 0;
         	a_gz = 0; 
				dabmode = 'off';
         end

	   	% rf excitation
	  		toppe.write2loop('tipdown.mod', 'RFphase', rfphs, 'RFamplitude', 1.0, 'textra', textra1);

		 	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
   		toppe.write2loop('readout.mod', 'DAQphase', rfphs, ...
				'slice', max(iz,1), 'echo', iim, 'view', iy, ...
				'dabmode', dabmode, 'Gamplitude', [1 a_gy a_gz]', 'textra', textra2);

		   % update rf phase (RF spoiling)
			rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
      end
	end
end
fprintf('\n');
toppe.write2loop('finish');

%% Save parameters and create tar file. We can reuse modules.txt from above.
save params seq seqb0
system('mkdir -p tar');
cd tar
system('tar cf scan,b0.tar ../main.m ../params.mat ../*.mod ../modules.txt ../scanloop.txt');
cd ..

fprintf('Scan time for 3D sequence: %.2f min\n', toppe.getscantime/60);

return;

