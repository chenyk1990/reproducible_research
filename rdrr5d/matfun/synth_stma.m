function synth_stma(stma1,stma2)
% Author      : Oboue et al. 2020
%               Zhejiang University
% 
% Date        : January, 2021
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) Oboue et al. 2020
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


%% load data
load ob_synth5d.mat
d=synth5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);
%
%% simultaneous denoising and reconstruction
randn('state',201314);
var=0.3;
dn=d+var*randn(size(d));

%% decimate
%[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.3;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0=dn.*mask;

%% stma 1st case 
K1=80; e1=0.999; T1=2; 
flow=0;fhigh=90;dt=0.004;Niter=15;mode=1;verb=1;iflb=0; 
a=(Niter-(1:Niter))/(Niter-1);
d1=denoising_reconstruction_stma5d(d0,mask,flow,fhigh,dt,Niter,eps,verb,mode,a,K1,e1,T1);

%% stma 2nd case 
K2=15; e2=0.999; T2=17; 
d2=denoising_reconstruction_stma5d(d0,mask,flow,fhigh,dt,Niter,eps,verb,mode,a,K2,e2,T2);

SNR(d, dn)
SNR(d, d0)
SNR(d, d1)
SNR(d, d2)

%% save data
% fid1=fopen('synth_clean5d.bin','w');
% fwrite(fid1,d,'float');
% fid2=fopen('synth_noisy5d.bin','w');
% fwrite(fid2,dn,'float');
% fid3=fopen('synth_obs5d.bin','w');
% fwrite(fid3,d0,'float');
% fid4=fopen('synth_rr5d.bin','w');
% fwrite(fid4,d1,'float');
% fid5=fopen('synth_irr5d.bin','w');
% fwrite(fid5,d2,'float');
%% from Matlab to Madagascar

%rsf_create(clean,size(d)');
%rsf_write(d,clean);

%rsf_create(noisy,size(dn)');
%rsf_write(dn,noisy);

%rsf_create(obs,size(d0)');
%rsf_write(d0,obs);

rsf_create(stma1,size(d1)');
rsf_write(d1,stma1);

rsf_create(stma2,size(d2)');
rsf_write(d2,stma2);






