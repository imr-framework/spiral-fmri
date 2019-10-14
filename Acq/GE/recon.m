  function imsos = recon(pfile, ksp, frames, params, sampWin, fovshift)
% function imsos = recon(pfile, ksp, frames, params, sampWin, fovshift)
%
% Combine data from 'frames' and reconstruct. Assumes fully-sampled stack-of-spirals data.
%
% Example:
%  >> p = getparams();
%  >> main(p);         % generates ksp.mat (and tar/scan,fmri.tar)
%  >> load ksp;
%  >> ims = recon('P12345.7', ksp, 13:15, params, 3:2770, [0 0 4]);

seq = params;

if ~exist('fovshift', 'var')
	fovshift = [0 0 0];
end

if length(frames) < seq.fmri.nLeafs
	warning('A fully sampled dataset requires at least nLeafs frames.');
end

% acquisition parameters
nLeafs  = seq.fmri.nLeafs;
fov     = seq.fov;        % cm
n       = seq.n;          % in-plane matrix size

% get data
[dat, rdb_hdr] = toppe.utils.loadpfile(pfile);    % [ndat ncoil nslice necho nview] = [ndat ncoil nslice 1 nframes]
dat = flipdim(dat,1);
%dat = circshift(dat,0);
dat = double(dat); 

% pick out one time-frame
dat = dat(:,:,:,:,frames);    % [ndat ncoil nslice 1 nLeafs]

% permute (as required by reconSoS below)
dat = permute(dat, [1 5 3 4 2]);         % [ndat nLeafs nslice 1 ncoils]            

% crop
kx = squeeze(ksp.kx(sampWin,1,1,frames));
ky = squeeze(ksp.ky(sampWin,1,1,frames));
dat = dat(sampWin,:,:,:,:);

% reconstruct
% dat = [ndat nleafs nz ntp ncoils]
% kx  = [ndat leafs]
kx = double(kx);
ky = double(ky);
[imsos] = toppe.utils.spiral.reconSoS(dat, kx, ky, [fov fov], [n n]);

% shift FOV
%ims(:,:,:,ifr) = circshift(imtmp, fovshift);

% iamge is mirrored about center. TODO: fix, or a 'feature' of nufft?
imsos = flipdim(imsos,1);
imsos = flipdim(imsos,2);

% display
%im(imsos); colormap gray;

return;
