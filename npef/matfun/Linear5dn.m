function Linear5dn(clean,noisy,obs,nobs,mmask,rr,drr,rrn,drrn)
% Author      : Yangkang Chen
% Date        : Oct, 2021
% 
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2021 The University of Texas at Austin
%  Copyright (C) 2021 Yangkang Chen
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
% 
%  Reference:   Chen et al., GJI, 2016; Chen et al., GJI, 2020; Chen et
%  al., GEO, 2021

% clc;clear;close all;


%% load data
% addpath(genpath('~/chenyk.data/cyk_small_dataset/'));

load yc_synth5d.mat
d=data5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);
dt=0.004;
%% exploring the data
%1) ploting CMP gather
%figure;imagesc(reshape(d(:,:,:,1,1),100,10*10));

%2) ploting common offset gather
%figure;imagesc(reshape(d(:,5,5,:,:),100,10*10));

%% simultaneous denoising and reconstruction
randn('state',201314);
var=0.1;
% var=0;
dn=d+var*randn(size(d));

%% decimate
[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.3;
% ratio=0.5;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0c=d.*mask;
d0=dn.*mask;

% 
%% reconstruct
flow=1;fhigh=100;dt=0.004;N=3;Niter=10;mode=0;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=rr5d_lb_recon(d0c,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);


flow=1;fhigh=100;dt=0.004;N=3;NN=3;Niter=10;mode=0;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d2=drr5d_lb_recon(d0c,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);


%% recon and denoise
flow=1;fhigh=100;dt=0.004;N=3;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1n=rr5d_lb_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);


flow=1;fhigh=100;dt=0.004;N=3;NN=3;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d2n=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

yc_snr(d(:),dn(:))
yc_snr(d(:),d0c(:))
yc_snr(d(:),d0(:))
yc_snr(d(:),d1(:))
yc_snr(d(:),d2(:))
yc_snr(d(:),d1n(:))
yc_snr(d(:),d2n(:))


%% from Matlab to Madagascar
rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0c)');
rsf_write(d0c,obs);

rsf_create(nobs,size(d0)');
rsf_write(d0,nobs);

rsf_create(mmask,size(mask)');
rsf_write(mask,mmask);

rsf_create(rr,size(d1)');
rsf_write(d1,rr);

rsf_create(drr,size(d2)');
rsf_write(d2,drr);

rsf_create(rrn,size(d1n)');
rsf_write(d1n,rrn);

rsf_create(drrn,size(d2n)');
rsf_write(d2n,drrn);

