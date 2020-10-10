function Hyper(clean,noisy,obs,sgk,ksvd,ddtf,mssa,mssasgk,T,niter,nniter,K,ll1,ll2,ll3,ss1,ss2,ss3,perc)
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

%% create data
cmp=zeros(50,32,32);
[nz,nx,ny]=size(cmp);

x=[1:32];
y=[1:32];
[x,y]=meshgrid(x,y);
z=10*(sqrt(1+((x-16)./10).^2+((y-16)./10).^2)-1)+20;
% figure;surf(y);set(gca,'ydir','reverse');

for ix=1:nx
    for iy=1:ny
        cmp(round(z(ix,iy))+1,ix,iy)=1;
    end  
end

cmp=reshape(cmp,nz,nx*ny);
wav=yc_ricker(30,0.004,0.2);
for ix=1:nx*ny
cmp(:,ix)=conv(cmp(:,ix),wav,'same');
end
cmp=reshape(cmp,nz,nx,ny);
dc=yc_scale(cmp,3);
%figure;imagesc(reshape(dc,50,32*32));colormap(seis);

%% adding noise
randn('state',201314);
var=0.1;
dn=dc+var*randn(size(dc));

% decimate
[nt,nx,ny]=size(dc);
ratio=0.5;
% mask=yc_genmask(reshape(dc,nt,nx*ny),ratio,'c',201415);
mask=ones(size(dc));
size(mask)
mask(:,2:2:end,:)=zeros(size(mask(:,2:2:end,:)));
size(mask)
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;
%figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10)]);colormap(seis);

%% simultaneous denoising and reconstruction
% Linear initialization (not useful ?)
% d1=yc_dlrecon_init(d0,mask,'nearest');
% % d1=InpaintingInterp2(d0,mask,'nearest');
% figure;imagesc([dc,d0,d1]);colormap(seis);


%benchmark with FXYMSSA
flow=0;fhigh=125;dt=0.004;N=25;Niter=12;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d4=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
%figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d4(:,:,10)]);colormap(seis);

%% SGK
% %param=struct('T',10,'niter',10,'mode',1,'K',128);
%mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=6;Niter=5; obtains 6.94 dB

%param=struct('T',10,'niter',10,'mode',1,'K',128);
%mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=6;Niter=5; %
param=struct('T',T,'niter',nniter,'mode',1,'K',K);
mode=1;l1=ll1;l2=ll2;l3=ll3;s1=ss1;s2=ss2;s3=ss3;perc=perc;Niter=niter; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
param.d0=d4;%param=rmfield(param,'d0');
% param.init=1;
tic;
d1=yc_sgk_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
%without initialization:  dB
%with d4 initialization:  dB (5 iterations)
toc;
yc_snr(dc,d1,2)
%figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d1(:,:,10)]);colormap(seis);

%% KSVD
% param=struct('T',2,'niter',10,'mode',1,'K',64);
% mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=1;Niter=10; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
tic;
%d2=yc_ksvd_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
d2=d0;
toc;
%figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d2(:,:,10)]);colormap(seis);
yc_snr(dc,d2,2)

%% DDTF
% param=struct('T',2,'niter',10,'mode',1,'K',64);
% mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=1;Niter=10; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
tic;
%d3=yc_ddtf_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
d3=d0;
toc;
%figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d3(:,:,10)]);colormap(seis);
yc_snr(dc,d3,2)

%% using d4 as initial model
param=struct('T',10,'niter',10,'mode',1,'K',128);
mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=6;Niter=5; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
param.d0=d4;%param=rmfield(param,'d0');
% param.init=1;
tic;
%d5=yc_sgk_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
d5=d0;
%without initialization: 9.73 dB
%with d4 initialization: 12.54 dB (5 iterations)
toc;
yc_snr(dc,d5,2)

yc_snr(dc,dn,2) % 6.92 dB
yc_snr(dc,d0,2) % 2.21 dB
yc_snr(dc,d1,2) % 9.54 dB 	%62.16	s
yc_snr(dc,d2,2) % dB	    %294.45 s
yc_snr(dc,d3,2) % dB		%262.41 s
yc_snr(dc,d4,2) % dB
yc_snr(dc,d5,2) % dB


rsf_create(clean,size(dc)');
rsf_write(dc,clean);

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

rsf_create(mssa,size(d4)');
rsf_write(d4,mssa);

rsf_create(mssasgk,size(d5)');
rsf_write(d5,mssasgk);
























