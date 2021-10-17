function LDRR(dobs,dmask,dout)
% Author      : Yangkang Chen
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

%%from Madagascar to Matlab
% create memory
addpath(genpath('~/chenyk/matlibcyk'));

%fid=fopen('../real5d_smaller/real5d_250_10_10_21_10.bin','r');
%fidm=fopen('../real5d_smaller/real5d_mask_250_10_10_21_10.bin','r');

%d=fread(fid,[250,10*10*21*10],'float');
%dm=fread(fidm,[250,10*10*21*10],'float');

load yc_field.mat;
load yc_field_mask.mat;

[nt,n1,n2,n3,n4]=size(Data5D);
d0=Data5D;d0=yc_scale(d0(:),1);
d0=reshape(d0,nt,n1,n2,n3,n4);
mask=Data5Dmask;

d=d0;
dm=mask;


d=reshape(d,250,10,10,21,10);
dm=reshape(dm,250,10,10,21,10);

D=zeros(250,20,20,21,10);
Dm=zeros(250,20,20,21,10);
D(:,1:2:end,1:2:end,:,:)=d;
Dm(:,1:2:end,1:2:end,:,:)=dm;
% for i=1:20
% D(:,1:2:end,1:2:end,:,i)=d(:,:,:,:,i);
% Dm(:,1:2:end,1:2:end,:,i)=dm(:,:,:,:,i);
% end

flow=0;fhigh=100;dt=0.004;N=30;NN=3;Niter=10;mode=1;verb=1;iflb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr5d_lb_recon(D,Dm,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

%% from Matlab to Madagascar
rsf_create(dobs,size(D)');
rsf_write(D,dobs);

rsf_create(dmask,size(Dm)');
rsf_write(Dm,dmask);

rsf_create(dout,size(d1)');
rsf_write(d1,dout);
    

return

