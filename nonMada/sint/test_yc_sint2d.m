%% DEMO script of yc_sint2d.m
%
% References
% [1] Chen, Y., Chen, X., Wang, Y. and Zu, S., 2019. The interpolation of sparse geophysical data. Surveys in Geophysics, 40(1), pp.73-105.
% [2] Wang, H., Chen, Y., Saad, O.M., Chen, W., Oboué, Y.A.S.I., Yang, L., Fomel, S. and Chen, Y., 2022. A Matlab code package for 2D/3D local slope estimation and structural filtering. Geophysics, 87(3), pp.F1–F14.
% [3] Huang, G., Chen, X., Li, J., Saad, O.M., Fomel, S., Luo, C., Wang, H. and Chen, Y., 2021. The slope-attribute-regularized high-resolution prestack seismic inversion. Surveys in Geophysics, 42(3), pp.625-671.
% [4] Huang, G., Chen, X., Luo, C. and Chen, Y., 2020. Geological structure-guided initial model building for prestack AVO/AVA inversion. IEEE Transactions on Geoscience and Remote Sensing, 59(2), pp.1784-1793.

clc;clear;close all;

%% load data
fid=fopen('hyper_zero.bin');
hyper0=fread(fid,[501,256],'float');

%% create mask/sampling operator
hyperm=ones(size(hyper0));
hyperm(find(hyper0==0))=0;

%% Slope estimation, read the reference [2]
dip=yc_dip2dmask(yc_bandpass(hyper0,0.004,0,10),hyperm);

%% sparse data interpolation
hyper_recon=yc_sint2d(hyper0,hyperm,dip,20,10,1,0.01);
figure;imagesc([hyper0,hyper_recon,hyper0-hyper_recon]);




