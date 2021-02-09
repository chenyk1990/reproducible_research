function syn_optshrinkdamped_stma3d(clean,noisy,obs,opshrinkdamped,opshrinkdampedstma)
% Author      : Oboue et al. 2020
%               Zhejiang University
% 
% Date        : January, 2021
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) Oboue et al. 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) Yapo Abole Serge Innocent Oboue and Yangkang Chen
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
%  KEY REFERENCE
%  Oboue et al., 2021, Robust damped rank-reduction method for simultaneous denoising and reconstruction of 5-D seismic data, Geophysics, 86, V71–V89.
% 
%  OTHER REFERENCES
%  [1] Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  [2] Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  [3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [4] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [5] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [6] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [7] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

%% generate synthetic data
%% generate synthetic data

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
d=plane3d/max(max(max(plane3d)));
%% without noise
dn=d;

%% decimate
[nt,nx,ny]=size(d);
ratio=0.5;
mask=yc_genmask(reshape(d,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;
%% reconstruct
% flow=0;fhigh=125;dt=0.004;N=3;NN=2;Niter=10;mode=0;verb=1;
% d1=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode);
% d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode);
% figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9),d(:,:,9)-d1(:,:,9)]);
% figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9),d2(:,:,9)]);
%% simultaneous denoising and reconstruction
% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));
d0=dn.*mask;

flow=0;fhigh=125;dt=0.004;N=3;NN=2;Niter=15;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
% d1=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
% d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);
d1=fxydmssa_recon_optshrink(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);

K=2.10; 
e=0.991; 
ws=1;
d2=fxydmssa_recon_optshrink_stma(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a,K,e,ws);

% figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9),d(:,:,9)-d1(:,:,9)]);
%figure;imagesc([d(:,:,9),d0(:,:,9),d3(:,:,9)]);

SNR_d1=SNR(d, d1)
SNR_d2=SNR(d, d2)

%% from Matlab to Madagascar
rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(opshrinkdamped,size(d1)');
rsf_write(d1,opshrinkdamped);

rsf_create(opshrinkdampedstma,size(d2)');
rsf_write(d2,opshrinkdampedstma);


