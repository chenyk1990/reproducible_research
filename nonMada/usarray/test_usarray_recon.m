% Demo script for teleseismic denoising and reconstruction
% as introduced in Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
% 
% This takes about 10 minutes
% 
% Written by Yangkang Chen
% Feb, 2018
% Modified on Dec, 2020
% 
% REFERENCES
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497?V506.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.

clc;clear;close all;
load('200901181411_wfm.mat');
% d, dists(shot receiver distance/offset in degree), stla, stlo, t 

%% rm bad trace
inds=[18,41,70];
d(:,inds)=[];d=yc_scale(d);
stlo(inds)=[];
stla(inds)=[];
d0=d(:,105:433);
stlo0=stlo(105:433);
stla0=stla(105:433);

%% 3D processing/reconstruction
mla=[33,49];
mlo=[-116,-102];
%binning
[d3d,x1,y1,mask]=yc_bin3d(d0,stlo0,stla0,16,28,mlo(1),mla(1),mlo(2),mla(2));
[stlo1,stla1]=meshgrid(x1,y1);
figure;plot(stlo0,stla0,'bv');hold on;
plot(stlo1(:),stla1(:),'r*');
print(gcf,'-depsc','-r200','grids.eps');

figure;imagesc(squeeze(mask(1,:,:))');colorbar;set(gca,'YDir','normal');

figure;
subplot(1,2,1);imagesc(squeeze(d3d(:,13,:)));caxis([-0.05,0.05]);colormap(gray);
subplot(1,2,2);imagesc(squeeze(d3d(:,:,10)));caxis([-0.05,0.05]);colormap(gray);


%% test if mask is correct
tt=(d3d.*mask-d3d);
norm(tt(:)) %0->correct

ratio=size(find(mask==1))/size(mask(:));
fprintf('Sampling ratio is %g\n',ratio);

figure;
subplot(1,2,1);imagesc(squeeze(mask(:,13,:)));caxis([0,1]);colormap(jet);colorbar;
subplot(1,2,2);imagesc(squeeze(mask(:,:,10)));caxis([0,1]);colormap(jet);colorbar;


%% global processing
flow=0;fhigh=0.5;dt=1;N=8;Niter=10;mode=1;verb=1;eps=0.00001;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=fxymssa_recon(d3d,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);

%% local processing
[n1,n2,n3]=size(d3d);
param.dt=dt;
param.flow=flow;
param.fhigh=fhigh;
param.N=N;
param.niter=Niter;
param.a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing;
param.eps=0.00001;
param.mode=1;
param.verb=1;
n1win=2000;n2win=n2;n3win=n3;
r1=0.5;r2=0.5;r3=0.5;
% Main program goes here !
d2=win3d_mask(@localfxymssa_recon, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);

param.amode=2;
d3=win3d_mask(@localfxymssa_recon_auto, mask, param, d3d, n1win, n2win, n3win, r1, r2, r3);

ilon=12;
dtest=squeeze(d3d(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 0.5, 1.2],'color','w');
yc_wigbh(dtest,stla1,t(1:5400),0.005);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
annotate;
print(gcf,'-depsc','-r400','us_lon1.eps');  

dtest=squeeze(d1(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 0.5, 1.2],'color','w');
yc_wigbh(dtest,stla1,t(1:5400),0.005);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Global','Fontsize',20);
annotate;
print(gcf,'-depsc','-r400','us_lon2.eps');  

dtest=squeeze(d2(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 0.5, 1.2],'color','w');
yc_wigbh(dtest,stla1,t(1:5400),0.005);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Local','Fontsize',20);
annotate;
print(gcf,'-depsc','-r400','us_lon22.eps');  

% dtest=squeeze(d3(:,ilon,:));
% figure('units','normalized','Position',[0.2 0.4 0.5, 1.2],'color','w');
% yc_wigbh(dtest,stla1,t(1:5400),0.005);
% ylim([31.5,50.2]);
% ylabel('Latitude (^o)','Fontsize',20);
% xlabel('Time (s)','Fontsize',20);
% title('Local Auto','Fontsize',20);
% annotate;
% print(gcf,'-depsc','-r400','us_lon222.eps');  


%% zoomed comparison
dtest=squeeze(d3d(1088:2587,ilon,8:12));
figure('units','normalized','Position',[0.2 0.4 0.5, 0.35],'color','w');
yc_wigbh2(dtest,stla1(8:12,1),t(1088:2587),8);
ylim([37,39]);xlim([500,1400]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
print(gcf,'-depsc','-r400','us_lon1_z.eps');  

dtest1=squeeze(d1(1088:2587,ilon,8:12));
figure('units','normalized','Position',[0.2 0.4 0.5, 0.35],'color','w');
yc_wigbh2(dtest1,stla1(8:12,1),t(1088:2587),8);
ylim([37,39]);xlim([500,1400]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Global','Fontsize',20);
print(gcf,'-depsc','-r400','us_lon2_z.eps');  

dtest2=squeeze(d2(1088:2587,ilon,8:12));
figure('units','normalized','Position',[0.2 0.4 0.5, 0.35],'color','w');
yc_wigbh2(dtest2,stla1(8:12,1),t(1088:2587),8);
ylim([37,39]);xlim([500,1400]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Local','Fontsize',20);
print(gcf,'-depsc','-r400','us_lon22_z.eps');  

% dtest3=squeeze(d3(1088:2587,ilon,8:12));
% figure('units','normalized','Position',[0.2 0.4 0.5, 0.35],'color','w');
% yc_wigbh2(dtest3,stla1(8:12,1),t(1088:2587),8);
% % yc_wigbh([dtest,dtest1],1:10,t(1088:2587),0.1);
% ylim([37,39]);xlim([500,1400]);
% ylabel('Latitude (^o)','Fontsize',20);
% xlabel('Time (s)','Fontsize',20);
% title('Local Auto','Fontsize',20);
% print(gcf,'-depsc','-r400','us_lon222_z.eps');  
