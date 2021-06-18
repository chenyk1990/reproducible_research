function [ D1 ] = yc_meanf(D,nfw,ifb,axis)
%MEANFCYK: mean filter along first or second axis for 2D profile
%  IN   D:   	intput data 
%       nfw:    window size
%       ifb:    if use padded boundary (if not, zero will be padded)
%       axis:    temporal sampling interval
%      
%  OUT   D1:  	output data
% 
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%
%  Example: dsp/test_meanfcyk_repeat.m
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
% 
% Reference:
% Chen, Y., 2020, Automatic microseismic event picking via unsupervised
% machine learning, GJI, 1750–1764
%
% Related reference:
% Chen, Y., 2018, Fast waveform detection for microseismic imaging using unsupervised machine learning, GJI, 1185-1199.
% Zhang et al., 2020, Convolutional Neural Networks for Microseismic Waveform Classification and Arrival Picking, Geophysics, WA227–WA240.
% Chen, et al., 2019, Automatic Waveform Classification and Arrival Picking Based on Convolutional Neural Network, ESS, 1244-1261.
% Saad and Chen, 2020, Automatic waveform-based source-location imaging using deep learning extracted microseismic signals, Geophysics, KS171–KS183.
% Qu et al., 2020, Automatic high-resolution microseismic event detection via supervised machine learning, GJI, 1881–1895.
%
% For earthquake detection and picking
% Saad and Chen, 2021, CapsPhase: Capsule Neural Network for Seismic Phase Classification and Picking, TGRS.
% Saad et al., 2021, SCALODEEP: A Highly Generalized Deep Learning Framework for Real-time Earthquake Detection, JGR, e2020JB021473
% Saad and Chen, 2020, Earthquake Detection and P-wave Arrival Time Picking using Capsule Neural Network, TGRS.

if nargin==0
 error('Input data must be provided!');
end

if nargin==1
 nfw=7;
 ifb=1;
 axis=2;
end;

if nargin==2
 ifb=1;
 axis=2;    
end

% nfw should be odd
if mod(nfw,2)==0
    nfw=nfw+1;
end

if axis==2
   D=D.'; 
end

[n1,n2]=size(D);
nfw2=(nfw-1)/2;

if ifb==1
    D=[flipud(D(1:nfw2,:));D;flipud(D(n1-nfw2+1:n1,:))];
else
    D=[zeros(nfw2,n2);D;zeros(nfw2,n2)];    
end

% output data
D1=zeros(n1,n2);
for i2=1:n2
   for i1=1:n1
      D1(i1,i2)=mean(D(i1:i1+nfw-1,i2)); 
   end 
end
if axis==2
    D1=D1.';
end
return
