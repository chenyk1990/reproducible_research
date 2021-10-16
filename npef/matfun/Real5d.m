function Real5d(obs,mmask,drr,obs2,mmask2,mmask2t)
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
load yc_field.mat;
load yc_field_mask.mat;

[nt,n1,n2,n3,n4]=size(Data5D);
d0=Data5D;d0=yc_scale(d0(:),1);
d0=reshape(d0,nt,n1,n2,n3,n4);
mask=Data5Dmask;

%% simultaneous denoising and reconstruction
%% reconstruct
flow=1;fhigh=100;dt=0.004;N=20;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
NN=3;
d1=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

%% from Matlab to Madagascar
rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(mmask,size(mask)');
rsf_write(mask,mmask);

rsf_create(drr,size(d1)');
rsf_write(d1,drr);

D=zeros(nt,n1*2,n2*2,n3*2,n4*2);
mask=zeros(size(D));
D(:,1:2:end,1:2:end,1:2:end,1:2:end)=d1;
mask(:,1:2:end,1:2:end,1:2:end,1:2:end)=ones(size(d1));
norm(D(:)-D(:).*mask(:))
rsf_create(obs2,size(D)');
rsf_write(D,obs2);

rsf_create(mmask2,size(mask)');
rsf_write(mask,mmask2);

mask=zeros(size(D));
mask(:,:,:,1:2:end,1:2:end)=ones(size(mask(:,:,:,1:2:end,1:2:end)));
rsf_create(mmask2t,size(mask)');
rsf_write(mask,mmask2t);


