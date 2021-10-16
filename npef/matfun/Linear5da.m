function Linear5da(obs,mmask,drr)
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

%% load data
load yc_synth5d.mat
d=data5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);
dt=0.004;

D=zeros(nt,nhx*2,nhy*2,nx,ny);
mask=zeros(size(D));
D(:,1:2:end,1:2:end,:,:)=d;
mask(:,1:2:end,1:2:end,:,:)=ones(size(d));
mask(:,:,:,2:2:end,:)=zeros(size(mask(:,:,:,2:2:end,:)));
mask(:,:,:,:,2:2:end)=zeros(size(mask(:,:,:,:,2:2:end)));
D=D.*mask;

%% reconstruct
flow=1;fhigh=100;dt=0.004;N=3;NN=6;Niter=10;mode=0;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
%d1=drr5d_lb_recon(D,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
d1=D;

%% from Matlab to Madagascar
rsf_create(obs,size(D)');
rsf_write(D,obs);

rsf_create(mmask,size(mask)');
rsf_write(mask,mmask);

rsf_create(drr,size(d1)');
rsf_write(d1,drr);


