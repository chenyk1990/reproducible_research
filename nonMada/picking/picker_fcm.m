function [ O, U3 ] = picker_fcm(X,q,nsta, nlta)
% picker_fcm: First-arrival picking using fuzzy c-means clustering algorithm
% 
% p,sigma,stalta are used for clustering (power, std, stalta)
% 
% Input:
%   q:    analysis window length (odd)
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
   q=15;nsta=5;nlta=20; %nsta and nlta should be a little longer to be smoother
end

nw2=(q-1)/2;

[n1,n2]=size(X);

if n1==1
    X=X';
   n1=n2;n2=1; 
end

O=zeros(1,n2);
U3=zeros(2,n1-2*nw2,n2);

for i2=1:n2
    xn=X(:,i2);
p0=yc_meanf(xn.^2,q,0,1)*q;p=p0(nw2+1:end-nw2);%p=p((q-1)/2:end-(q-1)/2);

a0=yc_meanf(xn,q,0,1)*q/q;a=a0(nw2+1:end-nw2);%a=a0((q-1)/2:end-(q-1)/2);

sigma0=yc_meanf((xn-a0).^2,q,0,1)*q/q;sigma=sigma0(nw2+1:end-nw2);%sigma=sigma((q-1)/2:end-(q-1)/2);

kurt0 = yc_meanf((xn-a0).^4,q,0,1)*q/q./(sigma0.^2)-3;kurt=kurt0(nw2+1:end-nw2);

stalta=yc_stalta(xn,nsta,nlta);stalta=stalta(nw2+1:end-nw2);
%kurt=sqrt(yc_meanf((xn-a0).^2,q,0,1)*q/m);sigma=kurt((q-1)/2:end-(q-1)/2);
stalta=yc_meanf(stalta,2*q,0,1);
%% SVD based method
p=yc_scale(p,1);
a=yc_scale(a,1);
sigma=yc_scale(sigma,1);
kurt=yc_scale(kurt,1);
stalta=yc_scale(stalta,1);

xt=xn(nw2+1:end-nw2);
xt=yc_scale(xt,1);
% figure;
% subplot(5,1,1);plot(p);ylim([-1,1]);title('power');
% subplot(5,1,2);plot(a);ylim([-1,1]);title('Mean');
% subplot(5,1,3);plot(sigma);ylim([-1,1]);title('STD');
% subplot(5,1,4);plot(kurt);ylim([-1,1]);title('Kurtosis');
% subplot(5,1,5);plot(stalta);ylim([-1,1]);title('STALTA');

%% improved version from [a,p,stalta]
% e=[a,p,sigma,stalta];
e=[a,p,stalta];%also try [a,2*p,stalta]
% figure(100);plot(stalta);pause(0.1);

n_clusters = 2;
randn('state',201617);
%[center,U,obj_fcn] = fcm(e, n_clusters);
[center,U,obj_fcn] = yc_fcm(e, n_clusters);

% figure;subplot(2,1,1);plot(U(1,:));
% subplot(2,1,2);plot(U(2,:));

%% select onset
eps=0.5;
if sum(abs(U(1,1:2*q))) < sum(abs(U(2,1:3*q)))
n_onset=min(find(U(1,:)>eps))+nw2;
else
n_onset=min(find(U(2,:)>eps))+nw2;   
end

O(i2)=n_onset;
U3(:,:,i2)=U;
end

end
