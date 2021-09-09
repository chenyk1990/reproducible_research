%  DEMO script for calculating local orthogonalization and its application in random noise attenuation 
%  
%  Copyright (C) 2016 Yangkang Chen
%  
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%  
%  Reference:   1. Random noise attenuation using local signal-and-noise orthogonalization
%               Chen and Fomel, 2015, Geophysics
%               2. Ground-Roll Noise Attenuation Using a Simple and Effective Approach Based on 
%               Local Band-Limited Orthogonalization, Chen et al., 2015, IEEE Geoscience and Remote Sensing Letters
%               3. Iterative deblending with multiple constraints based on shaping regularization,
%               Chen, 2015, IEEE Geoscience and Remote Sensing Letters
%               4. Orthogonalized morphological reconstruction for weak signal detection in micro-seismic monitoring:
%               Methodology, Huang et al., 2018, GJI
%               5. Surface-related multiple leakage extraction using local primary-and-multiple 
%               orthogonalization, Zhang et al., 2020, Geophysics
%               6. Non-stationary local signal-and-noise orthogonalization, Chen et al.,
%               2021, Geophysics
%               7. Local primary-and-multiple orthogonalization for leaked internal multiple crosstalk estimation and attenuation on full-wavefield migrated images
%               Zhang, et al., 2021, Geophysics
clc;clear;close all;

%% generate synthetic data
%This synthetic data was used in Huang et al., 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
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
d0=shot;

%% add noise
[n1,n2,n3]=size(d0);
randn('state',201415);
n=0.1*randn(n1,n2,n3);
dn=d0+n;%adding noise

%% Noise attenuation
d1=fxydmssa(dn,0,120,0.004,3,1,0);	%DMSSA (when damping factor =1, there are heavy damages)
noi1=dn-d1;

%% prepare paramters for ortho
rect=[10,10,10];
eps=0;
niter=20;
verb=1;

%% calculate local orthogonalization
[d2,noi2,low]=localortho(d1,noi1,rect,niter,eps,verb);

%% calculate local similarity
simi1=localsimi(d1,noi1,[5,5,5],niter,eps,verb);
simi2=localsimi(d2,noi2,[5,5,5],niter,eps,verb);

%% plot results
figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
subplot(1,2,1);imagesc([reshape(dn,n1,n2*n3);reshape(d1,n1,n2*n3);reshape(noi1,n1,n2*n3)]);caxis([-0.6,0.6]);colormap(jet);title('Initial denoising');
subplot(1,2,2);imagesc([reshape(dn,n1,n2*n3);reshape(d2,n1,n2*n3);reshape(noi2,n1,n2*n3)]);caxis([-0.6,0.6]);colormap(jet);title('Local orthogonalization');

figure;
subplot(2,1,1);imagesc(reshape(simi1,n1,n2*n3));caxis([0,1]);colormap(jet);title('Local similarity: Initial denoising');
subplot(2,1,2);imagesc(reshape(simi2,n1,n2*n3));caxis([0,1]);colormap(jet);title('Local similarity: Local orthogonalization');

fprintf('SNR of initial denoising is %g\n',str_snr(d0,d1,2));
fprintf('SNR of local orthogonalization is %g\n',str_snr(d0,d2,2));




