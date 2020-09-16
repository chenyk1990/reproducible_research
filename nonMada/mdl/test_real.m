clc;clear;close all;
%Demo for microseismic data denoising using 1D dictionary learning
%Prepared By Yangkang Chen and Hang Wang
%Dec, 09, 2018
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

load mdl_real.mat

%% patch size l1*l2
l1=16;l2=1;

c1=16;c2=32;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
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
% DCT=kron(dct,dct);%2D DCT dictionary (64,256)

%% plot the first 64 atoms
figure;
for ia=1:16
    subplot(4,4,ia);plot(dct(:,ia));
end


%% decompose the image into patches:
X=yc_patch(d,1,16,1,8,1);



%% OMP using DCT
nd=size(X,2);
K=5;
for i2=1:nd
    G(:,i2)=yc_omp0(D,X(:,i2),K);
end

%further constrain it to be sparser
G=yc_pthresh(G,'ph',1);
X2=D*G;

[n1,n2]=size(d);
d2=yc_patch_inv(X2,1,n1,n2,16,1,8,1);


figure('units','normalized');
imagesc(G);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_G.eps');


%% K-SVD
%This will be subsituted by independent subroutines
addpath(genpath('~/chenyk.data2/various/packages/toolbox'));

params.data = X;
params.Tdata = 5;%sparsity level
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
d1=yc_patch_inv(X1,1,n1,n2,16,1,8,1);
figure;%imagesc([d,d1,d-d1]);
subplot(3,1,1);plot(d);ylim([-0.4,0.4]);
subplot(3,1,2);plot(d1);ylim([-0.4,0.4]);
subplot(3,1,3);plot(d-d1);ylim([-0.4,0.4]);


%DWT
% d3 = wden(d,'rigrsure','s','sln',4,'db2');
% tmp=dwt(d,'db1','mode','sym');

%BP
d3=yc_bp(d,1,0.001,0.01,0.2,0.21);

%EMD
d4=yc_micro_emd(d);


dt=2;t=[0:n1-1]*dt;
figure('units','normalized','Position',[0.2 0.4 0.4, 0.8]);
subplot(5,1,1);plot(t,d,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(5,1,2);plot(t,d4,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);

subplot(5,1,3);plot(t,d3,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);

subplot(5,1,4);plot(t,d2,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(d)','Fontsize',16);

subplot(5,1,5);plot(t,d1,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
title('(e)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_dn.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.7]);
subplot(4,1,1);plot(t,d-d4,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(4,1,2);plot(t,d-d3,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);

subplot(4,1,3);plot(t,d-d2,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);

subplot(4,1,4);plot(t,d-d1,'b','linewidth',2);ylim([-0.5,0.5]);xlim(dt*[0,n1-1]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(d)','Fontsize',16);

xlabel('Time (ms)','Fontsize',16);
annotation('textarrow',[0.5,0.56],[0.45,0.41],'color','r','linewidth',2,'Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_n.eps');


figure('units','normalized');
imagesc(Gksvd0);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_Gksvd0.eps');

figure('units','normalized');
imagesc(Gksvd);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_Gksvd.eps');

figure('units','normalized');
imagesc(D);colormap(jet);caxis([-0.5,0.5]);%colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Atom NO','Fontsize',16);
title('Dictionary Atoms','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_Dct.eps');

figure('units','normalized');
imagesc(Dksvd);colormap(jet);caxis([-0.5,0.5]);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Atom NO','Fontsize',16);
title('Dictionary Atoms','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_Dksvd.eps');

figure('units','normalized');
imagesc(X);colormap(jet);%caxis([-0.5,0.5]);%colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_X.eps');

figure('units','normalized');
imagesc(X1);colormap(jet);%caxis([-0.5,0.5]);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','real_X1.eps');



%% plot the first 64 atoms
figure('units','normalized','Position',[0.2 0.4 0.7 0.8]);
for ia=1:30
    subplot(6,5,ia);plot([1:l1],D(:,ia),'b','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
    %     if ismember(ia,[1,6,11,16,21])
    %        subplot(4,4,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
    %     end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    ytickformat('%.1f');
    if ismember(ia,[1,6,11,16,21,26])
        ylabel('Amplitude','Fontsize',16);
    else
        set(gca,'yticklabel',[]);
        
    end
    
    if ismember(ia,[26,27,28,29,30])
        xlabel('Sample NO','Fontsize',16);
    else
        set(gca,'xticklabel',[]);
    end
    
    xlim([1,l1]);
end
print(gcf,'-depsc','-r300','real_atoms.eps');




figure('units','normalized','Position',[0.2 0.4 0.7 0.8]);
for ia=1:30
    subplot(6,5,ia);plot([1:l1],Dksvd(:,ia),'b','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
    %     if ismember(ia,[1,6,11,16,21])
    %        subplot(4,4,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
    %     end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    ytickformat('%.1f');
    if ismember(ia,[1,6,11,16,21,26])
        ylabel('Amplitude','Fontsize',16);
    else
        set(gca,'yticklabel',[]);
        
    end
    
    if ismember(ia,[26,27,28,29,30])
        xlabel('Sample NO','Fontsize',16);
    else
        set(gca,'xticklabel',[]);
    end
    
    xlim([1,l1]);
end
print(gcf,'-depsc','-r300','real_atoms1.eps');










