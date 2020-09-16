clc;clear;close all;
%% DEMO script for the Multitask dictionary learning algorithm
% Written by Yangkang Chen
% March, 2020
% 
% Key Reference:
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
%
% Some subroutines introduced in:
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717?1727.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2020, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 222, 1846?1863. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51?KS61.
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.

%% load data
load amir.mat;
dn=NoisyData;
dc=Data;
figure;imagesc([dc(:,:,15),dn(:,:,15),dn(:,:,15)-dc(:,:,15)]);

%% mapping to range [0,1] for non-negative MF
[dn1,minv,maxv]=Norm(dn);
[n1,n2,n3]=size(dn);

%% parameters
sigma=est_noise(dn1);
l1=45;l2=3;o1=1;o2=1;
r=round(1.15*l1*l2);
lambda=1*sigma*sqrt(2*log(r));
opts.maxit = 3000;
opts.tol = 1e-3;
opts.verbose = true;
% X=yc_patch(dn(:,:,1),1,l1,l2,o1,o2 );size(X) %[135,15048];

%% patching
patches=cell(n3,1);
for i3 = 1:n3
    patches{i3}=yc_patch(dn1(:,:,i3),1,l1,l2,o1,o2 );
%     size(patches{k})
end

%% Multi-task DL
[A, S] = MTSNMF(patches, r, lambda, opts);
% [A, S] = MTSNMF_MU(patches, r, LA, lambda, opts);

%% unpatching
d1=zeros(n1,n2,n3);
for i3 = 1:n3
    d1(:, :, i3)= yc_patch_inv( A{i3}*S,1,n1,n2,l1,l2,o1,o2 );
end

%% inverse normalization
d1=iNorm(d1,minv,maxv);

%% SNR
yc_snr(dc,dn,2)
yc_snr(dc,d1,2)

%% figure
figure;imagesc([dn(:,:,15),d1(:,:,15),dn(:,:,15)-d1(:,:,15)]);colormap(seis);

figure;
subplot(1,2,2);yc_wigb(dn(:,:,15)-d1(:,:,15));
subplot(1,2,1);yc_wigb(d1(:,:,15))

