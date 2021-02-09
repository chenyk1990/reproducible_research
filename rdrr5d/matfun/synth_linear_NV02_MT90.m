function synth_linear_NV02_MT90(clean,noisy,obs,drr,rdrr)
% Author      : Oboue et al. 2020
%               Zhejiang University
% 
% Date        : January, 2021
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) Oboue et al. 2020
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
%
%  KEY REFERENCE
%  Oboue et al., 2021, Robust damped rank-reduction method for simultaneous denoising and reconstruction of 5-D seismic data, Geophysics, 86, V71–V89.
% 
%  OTHER REFERENCES
%  [1] Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  [2] Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  [3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [4] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [5] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [6] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [7] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

%% load data
load ob_synth5d.mat
d=synth5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);

%% simultaneous denoising and reconstruction
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));

%% decimate
[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.1;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0=dn.*mask;

% DRR
flow=0;fhigh=100;dt=0.004;N=6;NN=2;Niter=10;mode=1;verb=1;iflb=0; 
a=(Niter-(1:Niter))/(Niter-1);
d1=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

% RDRR
K=2.5; 
e=0.891; 
T=1;
d2=rdrr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a,K,e,T);

snr_noisy=SNR(d, dn)
snr_obs=SNR(d, d0)
snr_drr=SNR(d, d1)
snr_rdrr=SNR(d, d2)

%% Estimate the rms error 
% Extract common-midpoint gather 
% Clean data
% CMP clean data
p = 1;
z = 1;
cmp_clean= d(:,:,:,p,z);
[nt,nx,ny]=size(cmp_clean);

% DRR data
% CMP DRR data
p = 1;
z = 1;
cmp_drr= d1(:,:,:,p,z);
[nt,nx,ny]=size(cmp_drr);

% DRR data
% CMP RDRR data
p = 1;
z = 1;
cmp_rdrr= d2(:,:,:,p,z);
[nt,nx,ny]=size(cmp_rdrr);

% rms error
rmse_cmp_drr = calc_rmse(cmp_clean,cmp_drr)
rmse_cmp_rdrr = calc_rmse(cmp_clean,cmp_rdrr)

% Extract common-offset gather 
% Clean data
% Offset clean data

x = 1;
t = 1;
offset_clean= squeeze(d(:,x,t,:,:));

% DRR data
% Offset DRR data

x = 1;
t = 1;
offset_drr= squeeze(d1(:,x,t,:,:));
[nt,nx,ny]=size(offset_drr);

% RDRR data
% Offset RDRR data

x = 1;
t = 1;
offset_rdrr= squeeze(d2(:,x,t,:,:));
[nt,nx,ny]=size(offset_rdrr);

% rms error
rmse_offset_drr = calc_rmse(offset_clean,offset_drr)
rmse_offset_rdrr = calc_rmse(offset_clean,offset_rdrr)


%% save data
% fid1=fopen('synth_clean5d.bin','w');
% fwrite(fid1,d,'float');
% fid2=fopen('synth_noisy5d.bin','w');
% fwrite(fid2,dn,'float');
% fid3=fopen('synth_obs5d.bin','w');
% fwrite(fid3,d0,'float');
% fid4=fopen('synth_rr5d.bin','w');
% fwrite(fid4,d1,'float');
% fid5=fopen('synth_irr5d.bin','w');
% fwrite(fid5,d2,'float');
%% from Matlab to Madagascar

rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(drr,size(d1)');
rsf_write(d1,drr);

rsf_create(rdrr,size(d2)');
rsf_write(d2,rdrr);






