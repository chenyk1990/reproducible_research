function Synth5d(clean,noisy,obs,sgk,ksvd,ddtf,T,niter,nniter,K,ll1,ll2,ll3,ss1,ss2,ss3,perc)
% Author      : Yangkang Chen
%               Zhejiang University
%         
% Date        : April, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Yangkang Chen
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

addpath(genpath('~/chenyk/matlibcyk'));
%clc;clear;close all;
% DEMO script for 3D denoising and reconstruction based on KSVD
% By Yangkang Chen
% March, 2020
a=4444

load('~/chenyk.data/cyk_small_dataset/synth5dcyk.mat');
d=data5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);
dt=0.004;
%% exploring the data
%1) ploting CMP gather
% figure;imagesc(reshape(d(:,:,:,1,1),100,10*10));

%2) ploting common offset gather
% figure;imagesc(reshape(d(:,5,5,:,:),100,10*10));

%% without noise
dn=d;

%% decimate
[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.3;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);

%% simultaneous denoising and reconstruction
% adding noise
randn('state',201314);
var=0.1;
dn=d+var*randn(size(d));
d0=dn.*mask;

%% SGK
param=struct('T',3,'niter',10,'mode',1,'K',512);
mode=1;l1=6;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;perc=0.5;Niter=5; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
% tic;
% d1=yc_sgk_recon5d(d0,mask,mode,[l1,l2,l3,l4,l5],[s1,s2,s3,s4,s5],perc,Niter,a,param);
% toc; %2256.75s for mode=1;l1=4;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;perc=1;Niter=10;


%% KSVD
% param=struct('T',2,'niter',10,'mode',1,'K',256);
% mode=1;l1=4;l2=4;l3=4;l4=4;l5=4;s1=3;s2=3;s3=3;s4=3;s5=3;perc=1;Niter=10; %
% a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
param=struct('T',3,'niter',10,'mode',1,'K',512);
mode=1;l1=6;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;perc=0.5;Niter=5; %
% tic
% d2=yc_ksvd_recon5d(d0,mask,mode,[l1,l2,l3,l4,l5],[s1,s2,s3,s4,s5],perc,Niter,a,param);
% toc %8766.97 s

%% DDTF
% param=struct('T',2,'niter',10,'mode',1,'K',256);
% mode=1;l1=4;l2=4;l3=4;l4=4;l5=4;s1=3;s2=3;s3=3;s4=3;s5=3;perc=1;Niter=10; %
% a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
param=struct('T',3,'niter',10,'mode',1,'K',512);
mode=1;l1=6;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;perc=0.5;Niter=5; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data

% tic
% d3=yc_ddtf_recon5d(d0,mask,mode,[l1,l2,l3,l4,l5],[s1,s2,s3,s4,s5],perc,Niter,a,param);
% toc % s


load syn5d_d1_l1_6_l4_s1_2_s2_perc_0_5_Niter_5.mat

yc_snr(d(:),d0(:),2) %-0.2079
yc_snr(d(:),dn(:),2) %-0.6590
yc_snr(d(:),d1(:),2) %6.58 (s1=2)
yc_snr(d(:),d2(:),2) %6.52
yc_snr(d(:),d3(:),2) %3.22 dB (s1=4)

rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(sgk,size(d1)');
rsf_write(d1,sgk);

rsf_create(ksvd,size(d2)');
rsf_write(d2,ksvd);

rsf_create(ddtf,size(d3)');
rsf_write(d3,ddtf);

























