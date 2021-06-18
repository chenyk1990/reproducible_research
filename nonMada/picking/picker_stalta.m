function [ O,R ] = picker_stalta(X,nsta, nlta)
% picker_fcm: First-arrival picking using STA/LTA algorithm
%
% Input:
%   nsta: short time period
%   nlta: long time period
%   STA/LTA
%   i>=nsta,
%   STA_i = (1/nsta) * sum_{j=i-nsta}^{i} a_j
%   LTA_i = (1/nlta) * sum_{j=i-nlta}^{i} a_j
%   i<nsta,STA_i=0,LTA_i=0
% Output:
%   On-set picker
%
% Copyright: Yangkang Chen, Jul, 2017
%
% Example: see test_micro_fcm.m
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


if nargin==1
    q=15;nsta=5;nlta=20;
end


[n1,n2]=size(X);

if n1==1
    X=X';
    n1=n2;n2=1;
end

O=zeros(1,n2);
R=zeros(n1,n2);

for i2=1:n2
    xn=X(:,i2);
    
    e=yc_stalta(xn,5,10);
    
    ref=yc_scale(e,1);
    
    n_onset=min(find(ref>0.98));
    
    O(i2)=n_onset;
    R(:,i2)=ref;
    
end

end
