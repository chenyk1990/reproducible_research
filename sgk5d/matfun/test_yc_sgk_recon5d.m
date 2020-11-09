clc;clear;close all;
% DEMO script for 5D denoising and reconstruction based on SGK
% By Hang Wang and Yangkang Chen
% April, 2020
% 
%  Reference:
%  Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE TGRS, doi: 10.1109/TGRS.2020.3030740
%  Chen, 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, GJI, 222, 1717??1727.

%% load data
load yc_synth5d.mat
d=data5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);
dt=0.004;

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
%or 
%2256.75s for mode=1;l1=4;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;perc=1;Niter=10;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
tic;
d1=yc_sgk_recon5d(d0,mask,mode,[l1,l2,l3,l4,l5],[s1,s2,s3,s4,s5],perc,Niter,a,param);
toc; 
figure;
subplot(3,1,1);imagesc(reshape(d(:,:,:,1,1),100,10*10));
subplot(3,1,2);imagesc(reshape(d0(:,:,:,1,1),100,10*10));
subplot(3,1,3);imagesc(reshape(d1(:,:,:,1,1),100,10*10));

figure;
subplot(3,1,1);imagesc(reshape(d(:,1,1,:,:),100,10*10));
subplot(3,1,2);imagesc(reshape(d0(:,1,1,:,:),100,10*10));
subplot(3,1,3);imagesc(reshape(d1(:,1,1,:,:),100,10*10));

yc_snr(d(:),d1(:),2) %5.4795


%% KSVD
param=struct('T',3,'niter',10,'mode',1,'K',512);
mode=1;l1=6;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;perc=0.5;Niter=5; %
tic
d2=yc_ksvd_recon5d(d0,mask,mode,[l1,l2,l3,l4,l5],[s1,s2,s3,s4,s5],perc,Niter,a,param);
toc %8766.97 s
figure;
subplot(3,1,1);imagesc(reshape(d(:,:,:,1,1),100,10*10));
subplot(3,1,2);imagesc(reshape(d0(:,:,:,1,1),100,10*10));
subplot(3,1,3);imagesc(reshape(d2(:,:,:,1,1),100,10*10));
yc_snr(d(:),d2(:),2) %4.84


%% DDTF
param=struct('T',3,'niter',10,'mode',1,'K',512);
mode=1;l1=6;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;perc=0.5;Niter=5; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data

tic
d3=yc_ddtf_recon5d(d0,mask,mode,[l1,l2,l3,l4,l5],[s1,s2,s3,s4,s5],perc,Niter,a,param);
toc % s
figure;
subplot(3,1,1);imagesc(reshape(d(:,:,:,1,1),100,10*10));
subplot(3,1,2);imagesc(reshape(d0(:,:,:,1,1),100,10*10));
subplot(3,1,3);imagesc(reshape(d3(:,:,:,1,1),100,10*10));


yc_snr(d(:),d0(:),2) 
yc_snr(d(:),dn(:),2)
yc_snr(d(:),d1(:),2) 
yc_snr(d(:),d2(:),2) 
yc_snr(d(:),d3(:),2) 













