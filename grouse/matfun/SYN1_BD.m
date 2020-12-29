function SYN1_BD(clean,noisy,dssa,dsst)
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



d=zeros(401,200);
dn=zeros(401,200);

%dn=zeros(512,50);
rsf_read(d,clean);
rsf_read(dn,noisy);
yc_snr(d,dn,2)

%% denoise
flow=0;fhigh=125;dt=0.004;N=3;verb=0;
tic
d1=fxymssa_grouse(dn,flow,fhigh,dt,N,1,verb);
toc
% figure;imagesc([d,dn,d1,dn-d1]);
% colormap(seis);
yc_snr(d,d1,2)
%about 12s
%snr:4.49

tic
d2=fxymssa_grouse(dn,flow,fhigh,dt,N,2,verb);
% [ d22,n22,ratio ] = localortho(d2,dn-d2,[20,20,1],50,0.01,1);
toc
% figure;imagesc([d,dn,d2,dn-d2]);colormap(seis);
% % figure;imagesc([d,dn,d22,dn-d22]);
% colormap(seis);
yc_snr(d,d2,2)
%about 6s
%about 6.46

tic
d3=fxymssa_grouse(dn,flow,fhigh,dt,N,3,verb);
toc
% figure;imagesc([d,dn,d3,dn-d3]);
% colormap(seis);
yc_snr(d,d3,2)
%about 1.5s
%SNR: 5.28 (error too much)


tic
d4=yc_fkt(dn,'ps',20);
toc
% figure;imagesc([d,dn,d4,dn-d4]);
% colormap(seis);
yc_snr(d,d4,2)
%about 0.05s
%SNR: 4.3

tic
d5=fxdecon(dn,flow,fhigh,dt,N,2,verb);
toc
% figure;imagesc([d,dn,d5,dn-d5]);
% colormap(seis);
yc_snr(d,d5,2)

% fid=fopen('../data/syn1_sp_ssa.bin','w');
% fwrite(fid,d1,'float');
% 
% fid=fopen('../data/syn1_sp_sst.bin','w');
% fwrite(fid,d2,'float');


rsf_create(dssa,size(d1)');
rsf_write(d1,dssa);

rsf_create(dsst,size(d2)');
rsf_write(d2,dsst);



