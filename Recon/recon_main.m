 clc
clear all
close all

%% Add dat2mat code 
[Filename,Pathname] = uigetfile('*.dat','Pick the raw data file');
twix_obj = mapVBVDVE(fullfile(Pathname,Filename));
image_obj = twix_obj{2}.image;%{2}.image;

sizeData = image_obj.sqzSize; %npoints, nchannels, kz*nframes

ksp_dat = squeeze(image_obj(:,:,:,1,1));

kz = 54;
nframes = 12;

ksp_dat = reshape(ksp_dat, size(ksp_dat,1), size(ksp_dat,2), kz, nframes);
ksp_dat = permute(ksp_dat,[1  3 4  2 ]);


%% Fully sampled recon
import toppe-master.*
fov = [24 24];           % cm
imsize = [72 72];
dx = fov(1)/imsize(1);   % voxel size (cm)

acqWin = 86:929;

nleafs = 3;  
frames = 1:(1+nleafs-1);    % combine these frames
frames = 7:9;
% Load data from Columbia
% 3 undersampled frames form one fully sampled frame
%datdir = '../../Data/ADNIphantom/';
%load([datdir '20190628_spiralPRESTO_ADNIphantom_dat.mat']);
dat = ksp_dat;

% Added by MMJ
dat = rfphs_compensate(dat);
% End of addition

dat = dat(acqWin,:,:,4,:);   %added by MMJ
%dat = dat(acqWin,:,frames,:);                       % [ndat, nz, nleafs, ncoils]
[ndat,nt,nz,ncoils] = size(dat);%added by MMJ
%[ndat,nz,nt,ncoils] = size(dat);
%dat = permute(dat, [1 3 2 4]);                      % [ndat, nleafs, nz, ncoils]
dat = reshape(dat, ndat, nleafs, nz, 1, ncoils);    % [ndat nleafs nz 1 ncoils]

% load kspace and convert to cycles/cm
ktrajdir = '../../Data/kspace-trajectory/';
load([ktrajdir 'full_kspace_trajectory.mat']);
ksp = ktraj_full(acqWin,:,frames,:);               % [ndat nz nleafs 3]
ksp = permute(ksp, [1 3 2 4]);                     % [ndat nleafs nz 3]
ksp = ksp(:,:,1,1:2)/max(abs(ksp(:)))/2;           % cycles/sample (approximately)
ksp = ksp/dx;                                      % cycles/cm

kx = ksp(:,:,1);
ky = ksp(:,:,2);

% dat: [ndat nleafs nz ntp ncoils]
% kx/ky: [ndat nleafs]
[imsos] = toppe.utils.spiral.reconSoS(dat, kx, ky, fov, imsize);
im(imsos)
