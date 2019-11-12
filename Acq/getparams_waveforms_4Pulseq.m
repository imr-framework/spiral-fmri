% Get acquisition parameters and RF/gradient waveforms for 3D
% stack-of-spirals fMRI scan
%
% Author: Marina Manso Jimeno
% Last Updated: 10/07/2019
%
%% Create a directory to store params and waveforms
mkdir assets
addpath assets

%% Get parameters
params = getparams; % Call the function in GE folder
% Add Siemens system limits for Pulseq (experimental values)
params.sysSiemens = struct('max_grad', 32, 'grad_unit','mT/m', 'max_slew',130, 'slew_unit','T/m/s', 'grad_raster_time',10e-6);
params.sysSiemens.gamma = 42576000;  % in Hz/T  %Determined from Pulseq - do not change
params.fatsat.freq_offset = -440; % Hz 

%% Get waveforms 
% fat saturation module
slThick = 1000;     % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse.
rf_fatsat = toppe.utils.rf.makeslr(params.fatsat.flip, slThick, params.fatsat.tbw, params.fatsat.dur, 0, 'ftype', 'ls', 'writeModFile', false, 'system', params.sysGE);

% PRESTO spoiler gradients
nCyclesSpoil = 0;        % should be balanced, since the PRESTO spoilers are played separately
% Will be placed on two axes for best (RF) spoiling.
% Derate slew rate a bit to keep PNS < 80% of stimulation threshold.
gspoil1 = toppe.utils.makecrusher(  params.fmri.nCyclesSpoil, params.dz, 0, 0.9*params.sys.maxSlew/sqrt(2), params.sys.maxGrad);
gspoil2 = toppe.utils.makecrusher(2*params.fmri.nCyclesSpoil, params.dz, 0, 0.7*params.sys.maxSlew/sqrt(2), params.sys.maxGrad);

% Tip-down (slab excitation) module
[rf,gex] = toppe.utils.rf.makeslr(params.rf.flip, params.rf.slabThick, params.rf.tbw, params.rf.dur, nCyclesSpoil,'ftype', params.rf.ftype, 'system', params.sys, 'writeModFile', false);
rf = [0*gspoil2(:); zeros(2,1); rf(:); 0*gspoil1(:)];
gx = [1*gspoil2(:); zeros(2,1); 0*gex(:);  -gspoil1(:)];
gy = [0*gspoil2(:); zeros(2,1); 0*gex(:); 0*gspoil1(:)];
gz = [1*gspoil2(:); zeros(2,1); gex(:);  -gspoil1(:)];
rf_tipdown = toppe.utils.makeGElength(rf);
gx_tipdown = toppe.utils.makeGElength(gx);
gy_tipdown = toppe.utils.makeGElength(gy);
gz_tipdown = toppe.utils.makeGElength(gz);

% Readout module (balanced stack-of spirals)
Router = params.fmri.nLeafs;
fovvd  = params.fov;
xresvd = params.n;
[g] = toppe.utils.spiral.genspivd2(fovvd, xresvd, Router, 0.99*params.sys.maxGrad, 0.99*params.sys.maxSlew*10, params.fmri.dsamp);
g = [0; 0; g(:)];           % add a couple of zeroes to make sure k=0 is sampled (?)

% make balanced and same length
gx = toppe.utils.makebalanced(real(g(:)), 'maxSlew', params.sys.maxSlew/sqrt(2), 'maxGrad', params.sys.maxGrad);
gy = toppe.utils.makebalanced(imag(g(:)), 'maxSlew', params.sys.maxSlew/sqrt(2), 'maxGrad', params.sys.maxGrad);
n = max(length(gx), length(gy));
gx = [gx; zeros(n-length(gx), 1)];
gy = [gy; zeros(n-length(gy), 1)];

% make it spiral-in
gx = flipud(gx);
gy = flipud(gy);

% add partition (kz) encoding trapezoids
gzamp = (1/params.sys.raster)/(params.sys.gamma*params.fovz);     % Gauss/cm
zarea = gzamp*params.nz*params.sys.raster;                 % Gauss/cm*sec
gpe = -toppe.utils.trapwave2(zarea/2, params.sys.maxGrad, params.sys.maxSlew, params.sys.raster*1e3);
gx1 = [0*gpe(:); zeros(2,1);   gx(:); 0*gpe(:)];
gy1 = [0*gpe(:); zeros(2,1);   gy(:); 0*gpe(:)];
gz1 = [  gpe(:); zeros(2,1); 0*gx(:);  -gpe(:)];

% write .mod file
gx_readout = toppe.utils.makeGElength(gx1);
gy_readout = toppe.utils.makeGElength(gy1);
gz_readout = toppe.utils.makeGElength(gz1);

%% Store params and waveforms
save('assets/params.mat','params')
save('assets/rf_fatsat.mat','rf_fatsat')
save('assets/rf_tipdown.mat','rf_tipdown')
save('assets/gx_tipdown.mat','gx_tipdown')
save('assets/gy_tipdown.mat','gy_tipdown')
save('assets/gz_tipdown.mat','gz_tipdown')
save('assets/gx_readout.mat','gx_readout')
save('assets/gy_readout.mat','gy_readout')
save('assets/gz_readout.mat','gz_readout')







