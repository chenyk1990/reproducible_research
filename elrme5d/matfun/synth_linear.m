function synth_linear(clean,noisy,obs,drr,elrme)
% Author      : Oboue and Chen, 2021
%               Zhejiang University
% 
% Date        : January, 2021
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2021 Zhejiang University
%  Copyright (C) Oboue and Chen, 2021
%  
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

%% simultaneous denoising and reconstruction
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));

%% decimate
[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.2;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0=dn.*mask;

% DRR
flow=0;fhigh=125;dt=0.004;N=3;NN=3;Niter=10;mode=1;verb=1; 
a=(Niter-(1:Niter))/(Niter-1);

tic
d1=drr5d_denoising_recon_synth(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);
toc

% ELRME
ws=1; e1=0.977; o1=3.25; b1=23; e0=0.971; o0=3.75; b0=10;
a=(Niter-(1:Niter))/(Niter-1);

tic
d2=elrme5d_denoising_recon_synth(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode, a, ws, e1, o1, b1, e0, o0, b0);
toc

%Normalized root square error (RSE)
rse1 = rse(d1,d)
rse2 = rse(d2,d)

% Signal-to-noise ratio (SNR)
snr_noisy=SNR(d, dn)
snr_obs=SNR(d, d0)
snr_drr=SNR(d, d1)
snr_rdrr=SNR(d, d2)

%% from Matlab to Madagascar

rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(drr,size(d1)');
rsf_write(d1,drr);

rsf_create(elrme,size(d2)');
rsf_write(d2,elrme);






