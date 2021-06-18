%% test of fuzzy c-means clustering algorithm in microseismic event picking
% by Yangkang Chen, 2017
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

%% FCM example
% data = load('fcmdata.dat');  % load some sample data
% n_clusters = 3;              % number of clusters
% [center,U,obj_fcn] = fcm(data, n_clusters);

clc;clear;close all;

%% parameters
q=25;nsta=5;nlta=20;

%% multi-channel (clean)
%!cp ~/chenyk.data2/various/cyksmall/mic_syn1.bin ./micro_data1.mat
load micro_data1.mat
d=yc_scale(d,2);


[n1,n2]=size(d);
[ O1,R ] = picker_stalta(d,nsta, nlta);
% figure;yc_wigb_p(d,100*ones(50,1),1,1:50,1:351,0.2,'ro');
figure;yc_wigb_p(d,O1);
  ylabel('Time(s)','Fontsize',16,'fontweight','bold');
  xlabel('Trace','Fontsize',16,'fontweight','bold');
  title('Noise','Fontsize',16,'fontweight','bold');
  set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold'); 
%   print(gcf,'-depsc','-r200','fig.eps');  

%% multi-channel (noisy)
randn('state',201517);
dn=d+0.2*randn(size(d));
yc_snr(d,dn)
figure;subplot(1,2,1);yc_wigb(d);
subplot(1,2,2);yc_wigb(dn);
[ O,U ] = picker_fcm(dn,q,nsta, nlta);
figure;subplot(1,2,1);yc_wigb_p(dn,O);

[ O1,R ] = picker_stalta(dn,nsta, nlta);
subplot(1,2,2);yc_wigb_p(dn,O1);






















