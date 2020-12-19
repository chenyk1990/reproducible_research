function Syn(clean,noisy,obs,mask,rr,drr)
% Author      : Yangkang Chen
% Date        : Dec, 2020
% 
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 The University of Texas at Austin
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
% 
%  Reference:   Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

%% load data

fid=fopen('syn_us_d_200_100_100.bin','r');
d=fread(fid,[200,100*100],'float');
d=reshape(d,200,100,100);

fid=fopen('syn_us_dn_200_100_100.bin','r');
dn=fread(fid,[200,100*100],'float');
dn=reshape(dn,200,100,100);

fid=fopen('syn_us_dm_200_100_100.bin','r');
dm=fread(fid,[200,100*100],'float');
dm=reshape(dm,200,100,100);

fid=fopen('syn_us_d0_200_100_100.bin','r');
d0=fread(fid,[200,100*100],'float');
d0=reshape(d0,200,100,100);


if norm(d0(:)-dn(:).*dm(:))~=0
	error('Data mismatch');
end

yc_snr(d,dn,2)
yc_snr(d,d0,2)

%% local RR
Niter=10;
param.dt=0.004;
param.flow=0;
param.fhigh=100;
param.N=3;
param.niter=Niter;
param.a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing;
param.eps=0.00001;
param.mode=1;
param.verb=1;

n1win=50;n2win=50;n3win=50;
r1=0.5;r2=0.5;r3=0.5;
%% Main program goes here !
d1=win3d_mask(@localfxymssa_recon, dm, param, d0, n1win, n2win, n3win, r1, r2, r3);
%or
%d1=win3d_mask(@localfxymssa_recon_auto, dm, param, d0, n1win, n2win, n3win, r1, r2, r3);

%% global RR
flow=0;fhigh=100;dt=0.004;N=9;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d2=fxymssa_recon(d0,dm,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
yc_snr(d,d1,2)	
yc_snr(d,d2,2)	


%% from Matlab to Madagascar
rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(mask,size(dm)');
rsf_write(dm,mask);

rsf_create(rr,size(d1)');
rsf_write(d1,rr);

rsf_create(drr,size(d2)');
rsf_write(d2,drr);








