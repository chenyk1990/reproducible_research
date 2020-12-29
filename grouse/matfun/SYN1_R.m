function SYN1_R(clean,noisy,dssa,dsst,dsnrs)
% Author      : Yangkang Chen
%               Zhejiang University
%         
% Date        : Dec, 2020
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




%% adding noise
% randn('state',201314);
% var=0.2;
% dn=d+var*randn(size(d));
d=zeros(401,200);
%dn=zeros(512,50);
%rsf_read(d,'../syn1/2Dsynth-clean.rsf');
fid=fopen('../matfun/clean.bin','r');
d=fread(fid,[401,200],'float');

randn('state',201314);
d=yc_scale(d,2);
var=0.3;
%rsf_read(dn1,'../syn1/syn.rsf');
n=yc_bp(var*randn(size(d)),0.004,0,5,30,40);
dn=d+n;

dn1=d+var*randn(size(d));

yc_snr(d,dn1,2)

%% denoise
flow=0;fhigh=125;dt=0.004;N=3;verb=0;
tic
d1=fxymssa_grouse(dn1,flow,fhigh/3,dt,N,1,verb);
toc
%figure;imagesc([d,dn1,d1,dn-d1]);
%colormap(seis);
yc_snr(d,d1,2)
%about 12s
%snr:4.49

tic
d2=fxymssa_grouse(d+0.2*randn(size(d)),flow,fhigh/3,dt,N,2,verb);
% [ d22,n22,ratio ] = localortho(d2,dn-d2,[20,20,1],50,0.01,1);
toc
%figure;imagesc([d,dn1,d2,dn-d2]);colormap(seis);
% figure;imagesc([d,dn1,d22,dn-d22]);
%colormap(seis);
yc_snr(d,d2,2)
%about 6s
%about 6.46

tic
d3=fxymssa_grouse(dn1,flow,fhigh,dt,N,3,verb);
toc
%figure;imagesc([d,dn1,d3,dn-d3]);
%colormap(seis);
yc_snr(d,d3,2)
%about 1.5s
%SNR: 5.28 (error too much)


tic
d4=yc_fkt(dn1,'ps',20);
toc
%figure;imagesc([d,dn1,d4,dn-d4]);
%colormap(seis);
yc_snr(d,d4,2)
%about 0.05s
%SNR: 4.3

tic
d5=fxdecon(dn1,flow,fhigh,dt,N,2,verb);
toc
%figure;imagesc([d,dn1,d5,dn-d5]);
%colormap(seis);
yc_snr(d,d5,2)

% dn=randn(600,40,40);
% flow=0;fhigh=125;dt=0.004;N=5;verb=0;
% tic
% d1=fxymssa_grouse(dn,flow,fhigh,dt,N,1,verb);
% toc
% 
% 
% tic
% d2=fxymssa_grouse(dn,flow,fhigh,dt,N,2,verb);
% toc

%% 
% var=1.0, -15.57,-5.07,-0.3854,0.413,-2.09,0.8074
% var=0.9, -14.66,-4.11,-0.13,0.99,-1.25,1.29
% var=0.8, -13.6354,-3.0024,0.167,1.66,-0.3132,1.86
% var=0.7, -12.47, -1.75, 0.61, 2.44, 0.74,2.54
% var=0.6, -11.136, -0.307,1.35, 3.39, 1.95,3.38
% var=0.5, -9.5530, 1.3299,2.4628, 4.55, 3.3648, 4.51
% var=0.4, -7.615,3.3323,4.18,6.608,5.05,   5.88
% var=0.3, -5.116,5.84,6.72,8.23,7.18,7.83
% var=0.2, -1.59,9.4412,10.9668,11.5567,10.0817,10.7162
% var=0.1, 4.42, 15.6564, 19.783, 17.6109, 14.6731,15.59

var=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
snr1=[4.42,-1.59,-5.116,-7.615,-9.553,-11.136,-12.47,-13.64,-14.66,-15.57];
snr2=[15.66,9.441,5.84,3.33,1.32,-0.307,-1.75,-3.00,-4.11,-5.07];
snr3=[19.783,10.96,6.72,4.18,2.46,1.35,0.61,0.17,-0.13,-0.39];
snrs=[var';snr1';snr2';snr3'];
%figure;plot(var,snr1,'k',var,snr2,'b',var,snr3,'r');
%legend('Raw','SSA','ST');

flow=0;fhigh=40;dt=0.004;N=3;verb=0;
tic
n1=fxymssa_grouse(n,flow,fhigh,dt,N,1,verb);
toc
%figure;imagesc([n,n1,n-n1]);
%colormap(seis);
%yc_snr(d,d1,2)
%about 12s
%snr:4.49

tic
n2=fxymssa_grouse(n,flow,fhigh,dt,N,2,verb);
toc
%figure;imagesc([n,n2,n-n2]);
%colormap(seis);
%yc_snr(d,d2,2)
%about 6s
%about 6.46

tic
n3=fxymssa_grouse(n,flow,fhigh,dt,N,3,verb);
toc
%figure;imagesc([n,n3,n-n3]);
%colormap(seis);
%yc_snr(d,d3,2)
%about 1.5s
%SNR: 5.28 (error too much)


tic
n4=yc_fkt(n,'ps',20);
toc
%figure;imagesc([n,n4,n-n4]);
%colormap(seis);
%yc_snr(d,d4,2)
%about 0.05s
%SNR: 4.3
% 
% 
% fid=fopen('../data/syn1_dn.bin','w');
% fwrite(fid,dn,'float');
% 
% fid=fopen('../data/syn1_c.bin','w');
% fwrite(fid,d,'float');
% 
% fid=fopen('../data/syn1_ssa.bin','w');
% fwrite(fid,d1,'float');
% 
% fid=fopen('../data/syn1_sst.bin','w');
% fwrite(fid,d2,'float');
% 
% fid=fopen('../data/syn1_snrs.bin','w');
% fwrite(fid,snrs,'float');
% %size:401*200, 4.698,4.605,0.405,0.02


rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(dssa,size(d1)');
rsf_write(d1,dssa);

rsf_create(dsst,size(d2)');
rsf_write(d2,dsst);

rsf_create(dsnrs,size(snrs)');
rsf_write(snrs,dsnrs);


