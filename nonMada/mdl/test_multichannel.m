% Demo for microseismic data denoising using 1D dictionary learning
% Prepared By Yangkang Chen and Hang Wang
% Dec, 09, 2018
%
% Key Reference
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51-KS61.
%
% For more details about dictionary learning and denoising, please refer to
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717?1727.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2020, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 222, 1846?1863. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51-KS61.
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.


clc;clear;close all;
addpath(genpath('~/chenyk.data2/various/packages/toolbox'));
%% prepare the data
load mdl_mc.mat

d=yc_scale(d,2);
% d=d(:,1:5:end);
[n1,n2]=size(d);

randn('state',20181920);
dn=d+0.3*randn(size(d));

figure;
subplot(1,2,1);yc_wigb(d);
subplot(1,2,2);yc_wigb(dn);


%% patch size l1*l2
l1=50;l2=1;

c1=50;c2=50;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
%% DCT dictionary (dctmtx will generates orthogonal transform)
dct=zeros(c1,c2);
for k=0:1:c2-1
    V=cos([0:1:c1-1]'*k*pi/c2);
    if k>0
        V=V-mean(V);
    end
    dct(:,k+1)=V/norm(V);
end
D=dct;

%% decompose the image into patches:
X=yc_patch(dn,1,l1,1,round(l1/2),1);

%% OMP using DCT
nd=size(X,2);
K=2;
for i2=1:nd
    G(:,i2)=yc_omp0(D,X(:,i2),K);
end

G=yc_pthresh(G,'ph',1);
X2=D*G;

[n1,n2]=size(d);
d2=yc_patch_inv(X2,1,n1,n2,l1,1,round(l1/2),1);

figure;yc_wigb(d2);

figure;imagesc(X);colormap(jet);
figure;imagesc(X2);colormap(jet);

%% K-SVD
addpath(genpath('~/chenyk/dblend_dl/toolbox/Dictionary'));
params.data = X;
params.Tdata = 2;%sparsity level
params.initdict=D;
params.dictsize = c2;
params.codemode = 'sparsity';
% params.codemode = 'error';
params.iternum = 30;
params.memusage = 'high';

rng(201819,'twister');
[Dksvd,Gksvd,err] = ksvd(params,'');

Gksvd0=Gksvd;
Gksvd=yc_pthresh(Gksvd0,'ph',1);
X1=Dksvd*Gksvd;

[n1,n2]=size(d);
d1=yc_patch_inv(X1,1,n1,n2,l1,1,round(l1/2),1);

figure;imagesc(X1);colormap(jet);

figure;
subplot(1,2,1);yc_wigb(d1);
subplot(1,2,2);yc_wigb(d2);

figure;yc_wigb([dn,d1,dn-d1]);
figure;yc_wigb([dn,d2,dn-d2]);

figure;
subplot(3,1,1);plot(d(:,1));ylim([-1.5,1.5]);
subplot(3,1,2);plot(dn(:,1));ylim([-1.5,1.5]);
subplot(3,1,3);plot(d1(:,1));ylim([-1.5,1.5]);

figure;imagesc(Dksvd);colormap(jet);
figure;imagesc(D);colormap(jet);

dt=2;t=[0:n1-1]*dt;
figure;
subplot(1,2,1);yc_wigb(d(:,1:4:end),1,1:13,t);
xlabel('Channel','Fontsize',16);
ylabel('Time (ms)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(1,2,2);yc_wigb(dn(:,1:4:end),1,1:13,t);
xlabel('Channel','Fontsize',16);
% ylabel('Time (ms)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);
print(gcf,'-depsc','-r300','syn2.eps');

%% single trace denoising comparison
ix=1;
dt=2;t=[0:n1-1]*dt;
figure('units','normalized','Position',[0.2 0.4 0.4, 0.4]);
subplot(2,1,1);plot(t,d(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(2,1,2);plot(t,dn(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
print(gcf,'-depsc','-r300','syn2_d.eps');

figure('units','normalized','Position',[0.2 0.4 0.4, 0.9]);
subplot(6,1,1);plot(t,d(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(6,1,2);plot(t,dn(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);

subplot(6,1,3);plot(t,d2(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);

subplot(6,1,4);plot(t,d1(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(d)','Fontsize',16);

subplot(6,1,5);plot(t,dn(:,ix)-d2(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(e)','Fontsize',16);

subplot(6,1,6);plot(t,dn(:,ix)-d1(:,ix),'b','linewidth',2);ylim([-1.5,1.5]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
title('(f)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
% annotation('textarrow',[0.6,0.66],[0.45,0.41],'color','r','linewidth',2,'Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','syn2_dn.eps');




%% other plots
figure('units','normalized');
imagesc(Gksvd0);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn2_Gksvd0.eps');

figure('units','normalized');
imagesc(G);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn2_G.eps');

figure('units','normalized');
imagesc(Gksvd);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn2_Gksvd.eps');

figure('units','normalized');
imagesc(D);colormap(jet);caxis([-0.2,0.2]);%colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Atom NO','Fontsize',16);
title('Dictionary Atoms','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn2_Dct.eps');

figure('units','normalized');
imagesc(Dksvd);colormap(jet);caxis([-0.4,0.4]);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Atom NO','Fontsize',16);
title('Dictionary Atoms','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn2_Dksvd.eps');

figure('units','normalized');
imagesc(X);colormap(jet);%caxis([-0.5,0.5]);%colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn2_X.eps');

figure('units','normalized');
imagesc(X1);colormap(jet);%caxis([-0.5,0.5]);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn2_X1.eps');




%% plot learned features
figure('units','normalized','Position',[0.2 0.4 0.8 2.0]);
for ia=1:50
    subplot(10,5,ia);plot([1:l1],D(:,ia),'b','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
    %     if ismember(ia,[1,2,3,4,5,13,14,28,37,43])
    %        subplot(10,5,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
    %     end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    ytickformat('%.1f');
    
    if ismember(ia,[1,6,11,16,21,26,31,36,41,46])
        ylabel('Amp','Fontsize',16);
            else
        set(gca,'yticklabel',[]);
        
        
    end
    
    if ismember(ia,[46,47,48,49,50])
        xlabel('Sample NO','Fontsize',16);
            else
        set(gca,'xticklabel',[]);
        
    end

    xlim([1,l1]);
end
print(gcf,'-depsc','-r300','syn2_atoms.eps');

figure('units','normalized','Position',[0.2 0.4 0.8 2.0]);
for ia=1:50
    subplot(10,5,ia);plot([1:l1],Dksvd(:,ia),'b','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
    if ismember(ia,[1,2,3,4,5,13,14,28,37,43])
        subplot(10,5,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    ytickformat('%.1f');
    
    if ismember(ia,[1,6,11,16,21,26,31,36,41,46])
        ylabel('Amp','Fontsize',16);
    else
        set(gca,'yticklabel',[]);
        
    end
    
    if ismember(ia,[46,47,48,49,50])
        xlabel('Sample NO','Fontsize',16);
    else
        set(gca,'xticklabel',[]);
        
    end
    
    xlim([1,l1]);
end
print(gcf,'-depsc','-r300','syn2_atoms1.eps');





