function synth(original,noisy,obs,hooi,drr,mrc)
% Author      : Oboue and Chen 2021
%               Zhejiang University
% 
% Date        : May, 2021
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2021 Zhejiang University
%  Copyright (C) Oboue et Chen 2021
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

%% load data
load yc_synth5d.mat
d=data5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);
dt=0.004;

%% simultaneous denoising and reconstruction
randn('state',201314);
var=0.1;
dn=d+var*randn(size(d));
%save d0.mat d0

%% decimate
[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.2;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0=dn.*mask;
%save d0.mat d0

%% HOOI 
flow=5;fhigh=100;dt=0.004;Niter=15;mode=1;verb=1;iflb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
rk=[3,3,3,3];

d1=hooi5d_lb_recon_synth(d0,mask,flow,fhigh,dt,rk,Niter,eps,verb,mode,2,a);
%save d1.mat d1

%% DRR
flow=5;fhigh=100;dt=0.004;r=3;NN=3;Niter=15;mode=1;verb=1;iflb=0; 
a=(Niter-(1:Niter))/(Niter-1);
d2=drr5d_lb_recon_synth(d0,mask,flow,fhigh,dt,r,NN,Niter,eps,verb,mode,iflb,a);
%save d2.mat d2

%% MRC
rk=[10,10,10,10];
d3=mrc5d_lb_recon_synth(d0,mask,flow,fhigh,dt,r,NN,rk,Niter,eps,verb,mode,iflb,a);
%save d3.mat d3

%% Signal-to-noise ratio (SNR)
snr_noisy=SNR(d, dn)
snr_obs=SNR(d, d0)
snr_HOOI=SNR(d, d1)
snr_drr=SNR(d, d2)
snr_mrc=SNR(d, d3)

%% from Matlab to Madagascar

rsf_create(original,size(d)');
rsf_write(d,original);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(hooi,size(d1)');
rsf_write(d1,hooi);

rsf_create(drr,size(d2)');
rsf_write(d2,drr);

rsf_create(mrc,size(d3)');
rsf_write(d3,mrc);






