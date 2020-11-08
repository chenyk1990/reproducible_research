function Field(clean,noisy,obs,sgk,ksvd,ddtf,mssa,mssasgk,T,niter,nniter,K,ll1,ll2,ll3,ss1,ss2,ss3,perc)
% Author      : Hang Wang and Yangkang Chen
%               Zhejiang University
%         
% Date        : April, 2020
%
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%
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
fid=fopen('./real3d.bin','r');
dc=fread(fid,[300,1000],'float');
dc=reshape(dc,300,100,10);

%% without noise
dn=dc;

%% decimate
[nt,nx,ny]=size(dc);
ratio=0.5;
mask=yc_genmask(reshape(dc,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;

%% simultaneous denoising and reconstruction
% adding noise
randn('state',201314);
var=0.1;
dn=dc+var*reshape(bpn(nt,nx*ny,0.004,0,1,35,36,201920),nt,nx,ny);
d0=dn.*mask;

%benchmark with FXYMSSA
flow=0;fhigh=125;dt=0.004;N=10;Niter=10;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d4=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);

%% SGK
param=struct('T',T,'niter',nniter,'mode',1,'K',K);
mode=1;l1=ll1;l2=ll2;l3=ll3;s1=ss1;s2=ss2;s3=ss3;perc=perc;Niter=niter; %
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing for noisy data
param.d0=d4;%param=rmfield(param,'d0');
tic;
d1=yc_sgk_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
toc;
yc_snr(dc,d1,2)
%figure;imagesc([dc(:,:,10),dn(:,:,10),d0(:,:,10),d1(:,:,10)]);colormap(seis);

%% KSVD
tic;
d2=yc_ksvd_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
toc;
yc_snr(dc,d2,2)

%% DDTF
tic;
d3=yc_ddtf_recon(d0,mask,mode,[l1,l2,l3],[s1,s2,s3],perc,Niter,a,param);
toc;
yc_snr(dc,d3,2)

%% using d4 as initial model
param=struct('T',10,'niter',10,'mode',1,'K',128);
mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=6;Niter=5; %
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
























