% Author      : Hang Wang and Yangkang Chen
%               Zhejiang University
%         
% Date        : April, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Hang Wang and Yangkang Chen
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
%  Reference:
%  Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE TGRS, doi: 10.1109/TGRS.2020.3030740
%  Chen, 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, GJI, 222, 1717??1727.

%% create data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
dc=yc_scale(plane3d,3);
dc=dc(51:225,:,:);
%figure;imagesc(reshape(dc,175,20*20));colormap(seis);

%% adding noise
randn('state',201314);
var=0.1;
dn=dc+var*randn(size(dc));

% decimate
[nt,nx,ny]=size(dc);
ratio=0.5;
mask=yc_genmask(reshape(dc,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;
figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10)]);colormap(seis);

%% SGK
param=struct('T',3,'niter',10,'mode',1,'K',64);
mode=1;l1=6;l2=4;l3=4;s1=2;s2=2;s3=2;perc=4;Niter=12; %
% param=struct('T',T,'niter',nniter,'mode',1,'K',K);
% mode=1;l1=ll1;l2=ll2;l3=ll3;s1=ss1;s2=ss2;s3=ss3;perc=perc;Niter=niter; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
tic;
d1=yc_sgk_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
toc;
yc_snr(dc,d1,2)
figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d1(:,:,10)]);colormap(seis);


%% KSVD
% param=struct('T',2,'niter',10,'mode',1,'K',64);
% mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=1;Niter=10; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
% a=ones(Niter,1);
tic;
d2=yc_ksvd_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
toc;
figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d2(:,:,10)]);colormap(seis);
yc_snr(dc,d2,2)

%% DDTF
% param=struct('T',2,'niter',10,'mode',1,'K',64);
% mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=1;Niter=10; %
tic;
d3=yc_ddtf_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
toc;
figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d3(:,:,10)]);colormap(seis);
yc_snr(dc,d3,2)

%benchmark with FXYMSSA
flow=0;fhigh=125;dt=0.004;N=3;Niter=12;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d4=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d4(:,:,10)]);colormap(seis);

%% using d4 as initial model
param=struct('T',3,'niter',10,'mode',1,'K',64);
mode=1;l1=6;l2=4;l3=4;s1=2;s2=2;s3=2;perc=5;Niter=5; %
param.d0=d4;%param=rmfield(param,'d0');
tic;
d5=yc_sgk_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
toc;
yc_snr(dc,d5,2)

yc_snr(dc,dn,2) 
yc_snr(dc,d0,2)
yc_snr(dc,d1,2) 
yc_snr(dc,d2,2) 
yc_snr(dc,d3,2) 
yc_snr(dc,d4,2) 
yc_snr(dc,d5,2) 


tic
[d11,D1,G1,DCT1]=yc_sgk_denoise(d1,mode,[l1,l2,l3],[s1,s2,s3],perc,param);
toc
tic
[d22,D2,G2,DCT2]=yc_ksvd_denoise(d2,mode,[l1,l2,l3],[s1,s2,s3],perc,param);
toc
tic
[d33,D3,G3,DCT3]=yc_ddtf_denoise(d3,mode,[l1,l2,l3],[s1,s2,s3],perc,param);
toc




















