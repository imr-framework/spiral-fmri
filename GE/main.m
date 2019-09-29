function main(seq)
% Create stack-of-spirals PRESTO fMRI sequence for TOPPE.
% Also create matched 3D spin-warp B0 sequence (also used for sensitivity map estimation)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First create stack-of-spirals fMRI sequence                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create balanced stack-of-spirals readout leaf (readout.mod). Isotropic resolution.
% design spiral
rmax = seq.n/(2*seq.fov);           % max k-space radius
[~,g] = toppe.utils.spiral.mintgrad.vds(0.99*seq.sys.maxSlew*1e3, 0.99*seq.sys.maxGrad, seq.sys.raster, seq.fmri.nLeafs, seq.fov, 0, 0, rmax);   % vds returns complex k, g
g = [0; 0; g(:)];           % add a couple of zeroes to make sure k=0 is sampled

% make balanced and same length
gx = toppe.utils.makebalanced(real(g(:)), 'maxSlew', seq.sys.maxSlew/sqrt(2), 'maxGrad', seq.sys.maxGrad);
gy = toppe.utils.makebalanced(imag(g(:)), 'maxSlew', seq.sys.maxSlew/sqrt(2), 'maxGrad', seq.sys.maxGrad);
n = max(length(gx), length(gy));
gx = [gx; zeros(n-length(gx), 1)];
gy = [gy; zeros(n-length(gy), 1)];

% make it spiral-in
gx = flipud(gx);
gy = flipud(gy);

% add partition (kz) encoding trapezoids
gzamp = (1/seq.sys.raster)/(seq.sys.gamma*seq.fovz);     % Gauss/cm
zarea = gzamp*seq.nz*seq.sys.raster;                 % Gauss/cm*sec
gpe = -toppe.utils.trapwave2(zarea/2, seq.sys.maxGrad, seq.sys.maxSlew, seq.sys.raster*1e3);
gx1 = [0*gpe(:); zeros(2,1);   gx(:); 0*gpe(:)];
gy1 = [0*gpe(:); zeros(2,1);   gy(:); 0*gpe(:)];
gz1 = [  gpe(:); zeros(2,1); 0*gx(:);  -gpe(:)];

% write to .mod file
gx1 = toppe.utils.makeGElength(gx1);
gy1 = toppe.utils.makeGElength(gy1);
gz1 = toppe.utils.makeGElength(gz1);
toppe.writemod('gx', gx1, 'gy', gy1, 'gz', gz1, 'ofname', 'readout.mod', ...
               'desc', 'stack-of-spirals readout module', 'system', seq.sysGE);


%% Create slab-selective excitation (tipdown.mod). Include (PRESTO) spoiler gradients.

nCyclesSpoil = 0;       % make it balanced initially -- gradient pre/rephasers will then be modified
[rf,gex] = toppe.utils.rf.makeslr(seq.rf.flip, seq.rf.slabThick, seq.rf.tbw, seq.rf.dur, nCyclesSpoil, ...
                       'ftype', seq.rf.ftype, 'system', seq.sys, 'writeModFile', false);

% Create spoiler (PRESTO) gradients.
% Will be placed on two axes for best (RF) spoiling.
gspoil1 = toppe.utils.makecrusher(seq.fmri.nCyclesSpoil, seq.dz, 0, seq.sys.maxSlew, seq.sys.maxGrad);
gspoil2 = toppe.utils.makecrusher(2*seq.fmri.nCyclesSpoil, seq.dz, 0, seq.sys.maxSlew, seq.sys.maxGrad);

% create tipdown.mod
rf = [0*gspoil2(:); zeros(2,1); rf(:); 0*gspoil1(:)];
gx = [1*gspoil2(:); zeros(2,1); 0*gex(:);  -gspoil1(:)];
gy = [0*gspoil2(:); zeros(2,1); 0*gex(:); 0*gspoil1(:)];
gz = [1*gspoil2(:); zeros(2,1); gex(:);  -gspoil1(:)];
rf = toppe.utils.makeGElength(rf);
gx = toppe.utils.makeGElength(gx);
gy = toppe.utils.makeGElength(gy);
gz = toppe.utils.makeGElength(gz);
toppe.writemod('rf', rf, 'gx', gx, 'gy', gy, 'gz', gz, 'ofname', 'tipdown.mod', ...
               'desc', 'RF slab excitation with PRESTO gradients', 'system', seq.sysGE);  % NB! Pass 'sysGE' here, not 'sys'!

return;

%% Create fat saturation pulse
% bw = 500 Hz. Frequency offset (-440 Hz) is set in scanloop.txt.
seq.fatsat.flip = 50;
slThick = 1000;     % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse
seq.fatsat.tbw = 1.5;
seq.fatsat.dur = 3;            % pulse duration (msec)
toppe.utils.rf.makeslr(seq.fatsat.flip, slThick, seq.fatsat.tbw, seq.fatsat.dur, 0, ...
                       'ftype', 'ls', 'ofname', 'fatsat.mod', 'system', seq.sysGE);

%% Create scanloop.txt
nz_samp = 30;                  % Sample this many kz points per time-frame (undersampling factor in kz is 54/30 = 1.8)
seqTR = 16.7e-3;               % sequence TR (sec). See toppe.getTRtime()
dur = 5*60;     % total duration of fMRI scan (sec)
trVol = nz_samp*seqTR;         % time to acquire one under-sampled image volume (sec)
if test
	seq.nt = 30;
else
	seq.nt = 2*round(dur/trVol/2);      % number of (undersampled) time-frames
end

% fully sampled kz sampling pattern
for ii = 1:seq.nz
	kzFull(ii) = ((ii-1+0.5)-seq.nz/2)/(seq.nz/2);    % scaling is (-1 1)
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
	kzU = kzFull(logical(zInd));
end

seq.nframes = seq.nref + seq.nt;

rfphs = 0;              % radians
rfphsLast = rfphs;
daqphs = 0;
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 150;  % For 30 kz platters per frame and rf_spoil_seed=150, ...
                      % we have mod(nz_samp*rf_spoil_seed,360)=180 which is what we need for improved RF spoiling...
                      % using our method in Magn Reson Med. 2016 Jun;75(6):2388-93. doi: 10.1002/mrm.25843.

% loop over (undersampled) time frames and fill in scanloop.txt,
% and write k-space values to file
if writeKspace
	[rf,gx,gy] = toppe.readmod('readout.mod');  % one spiral leaf
	[kx1,ky1] = toppe.utils.g2k([gx(:,1) gy(:,1)],1);
	k1 = complex(kx1,ky1);  % single spiral 'prototype'
	ndat = size(k1,1);
	necho = 1;
	ksp.kx = NaN*ones(ndat, seq.nz, necho, seq.nframes);   % to match the 'slice-echo-view' order of 'dat' array returned by toppe.utils.loadpfile
	ksp.ky = ksp.kx;
	ksp.kz = ksp.kx;
end
fprintf('Writing scanloop.txt for fMRI sequence\n');
toppe.write2loop('setup');
for iframe = (-ndisdaq+1):seq.nframes
	if ~mod(iframe,10)
		fprintf([repmat('\b',1,20) sprintf('%d of %d', iframe+ndisdaq, seq.nt+ndisdaq+seq.nref )]);
	end

	% Set kz sampling pattern for this frame.
	if iframe < (seq.nref+1)
		kz1 = kzFull;                % numel(kz) = seq.nz;
	else
		kz1 = kzU;                   % numel(kz) = nz_samp;
	end

	% set 'view' data storage index
	if iframe < 1 
		dabmode = 'off';
		view  = 1;
	else
		dabmode = 'on';
		view = iframe;
	end

	for iz = 1:numel(kz1)
   	% fat sat
   	% toppe.write2loop('fatsat.mod', 'RFphase', rfphs, 'Gamplitude', [0 0 0]');

   	% rf excitation (and PRESTO gradients)
   	toppe.write2loop('tipdown.mod', 'RFphase', rfphs, 'Gamplitude', [1 0 1]');

   	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
		slice = iz;
		view = max(iframe,1);
		if iframe < (seq.nref+1)
			ileaf = mod(mod(iframe,seq.nleafs) + iz,seq.nleafs) + 1;   % rotate leaf every frame and every kz platter
		else
			ileaf = mod(iz,seq.nleafs) + 1;                        % rotate leaf every kz platter (i.e., same undersampling pattern for every frame)
		end
		%ileaf = mod(mod(iframe,seq.nleafs) + iz,seq.nleafs) + 1;   % rotate leaf every frame and every kz platter
		phi = 2*pi*(ileaf-1)/seq.nleafs;                          % leaf rotation angle (radians)
		echo = 1;
   	toppe.write2loop('readout.mod', 'DAQphase', rfphsLast, 'slice', slice, 'echo', echo, ...
			'view', view, 'Gamplitude', [1 1 kz1(iz)]', 'dabmode', dabmode, 'rot', phi);

	   % update rf phase (RF spoiling)
		rfphsLast = rfphs;
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;

		% kspace info for this TR
		if strcmp(dabmode, "on") & writeKspace
			k1tmp = k1.*exp(1i*phi)*seq.fov/seq.n;         % convert from cycles/cm to cycles/sample
			ksp.kx(:,slice,echo,view) = real(k1tmp);
			ksp.ky(:,slice,echo,view) = imag(k1tmp);
			ksp.kz(:,slice,echo,view) = kz1(iz)/2;         % cycles/sample
		end

		ksp.kinfo(slice,echo,view).ileaf = ileaf;     % might come in handy
		ksp.kinfo(slice,echo,view).rot = phi;         % this too
	end
end
fprintf('\n');
toppe.write2loop('finish');

%% Save k-space trajectory and other sequence info 
%[rf,gx,gy,gz,desc,paramsint16,paramsfloat] = toppe.readmod('readout.mod');
%[ksp.kx ksp.ky] = toppe.utils.g2k([gx(:) gy(:)], seq.nleafs);
if writeKspace
	ksp.kx = single(ksp.kx);
	ksp.ky = single(ksp.ky);
	ksp.kz = single(ksp.kz);
	fprintf('Writing ksp.mat...');
	%save -v7.3 ksp ksp
	save ksp ksp
	fprintf('done\n');
end

%% create tar file
system('mkdir -p tar');
cd tar
system('tar cf scan,fmri.tar ../main.m ../params.mat ../ksp.mat ../*.mod ../modules.txt ../scanloop.txt');
cd ..

%% display sequence
%toppe.playseq(2, 'tpause', 0.1);

%% convert to Pulseq
if 0
% addpath ~/gitlab/toppe/pulseq/
cd tar
ge2seq('scan.tar');
seq = mr.Sequence();
seq.read('out.seq');
seq.plot_sg('TimeRange', [10 10.05]);
cd ..
end

fprintf('Scan time for fMRI sequence: %.2f min\n', toppe.getscantime/60);

%keyboard


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

