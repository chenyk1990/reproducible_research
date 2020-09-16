clc;clear;close all;
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



addpath(genpath('~/chenyk.data2/various/packages/toolbox'));
r=zeros(256,1);
r(100)=1;
r(180)=0.7;
% r(240)=0.3;
w=yc_ricker(20,0.004,0.1);
d=conv(r,w,'same');
d=yc_phase_correction(d,90);
figure;plot(d);

d0=d;
randn('state',201516);
d=d0+0.1*randn(size(d0));

%% patch size l1*l2
l1=16;l2=1;

c1=16;c2=30;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
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

% %% plot the first 64 atoms
% figure;
% for ia=1:16
%     subplot(4,4,ia);plot(dct(:,ia));
% end



%% decompose the image into patches:
X=yc_patch(d,1,16,1,8,1);

%% for ppt
figure;plot(70:150,d(70:150),'r','linewidth',4);
Dt=X(:,3:14);Dt(:,2:2:end)=zeros(size(Dt(:,2:2:end)));
figure;yc_wigbh(Dt);
figure('units','normalized','Position',[0.2 0.4 0.3, 0.8],'color','w');wigbh_tmp(X(:,9:14)+0.25*randn(16,6),1:6,1:16,1.5);

ylim([0,7]);


%% OMP using DCT
nd=size(X,2);
K=2;
for i2=1:nd
   G(:,i2)=yc_omp0(D,X(:,i2),K); 
end

%further constrain it to be sparser
G=yc_pthresh(G,'ph',0.4);
X2=D*G;

[n1,n2]=size(d);
d2=yc_patch_inv(X2,1,n1,n2,16,1,8,1);
yc_snr(d0,d2)   %6.63

figure;%imagesc([d,d1,d-d1]);
subplot(3,1,1);plot(d);ylim([-1,1]);
subplot(3,1,2);plot(d2);ylim([-1,1]);
subplot(3,1,3);plot(d-d2);ylim([-1,1]);

figure('units','normalized');
imagesc(G);colormap(jet);caxis([-1,2.0]);colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_G.eps');

%% K-SVD
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
Gksvd=yc_pthresh(Gksvd0,'ph',0.4);
X1=Dksvd*Gksvd;

[n1,n2]=size(d);
d1=yc_patch_inv(X1,1,n1,n2,16,1,8,1);
dt=2;t=[0:n1-1]*dt;

t1=170;t2=225;t3=330;t4=385;
figure('units','normalized','Position',[0.2 0.4 0.4, 0.4]);
subplot(2,1,1);plot(t,d0,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);

subplot(2,1,2);plot(t,d,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);
print(gcf,'-depsc','-r300','syn_d.eps');

%BP
d3=yc_bp(d,0.002,0,0.5,100,105);

%EMD
d4=yc_micro_emd(d);  

yc_snr(d0,d2)    %6.63
yc_snr(d0,d1)   %12.06
yc_snr(d0,d3)   %7.81
yc_snr(d0,d4)   %6.36

 figure('units','normalized','Position',[0.2 0.4 0.4, 0.6]);
subplot(4,1,1);plot(t,d4,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);
title('(a)','Fontsize',16);

subplot(4,1,2);plot(t,d3,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);
title('(b)','Fontsize',16);

subplot(4,1,3);plot(t,d2,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);
title('(c)','Fontsize',16);

subplot(4,1,4);plot(t,d1,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(d)','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);
print(gcf,'-depsc','-r300','syn_dn1.eps');

figure('units','normalized','Position',[0.2 0.4 0.4, 0.6]);
subplot(4,1,1);plot(t,d-d4,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);

subplot(4,1,2);plot(t,d-d3,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
title('(b)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);

 subplot(4,1,3);plot(t,d-d2,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);

subplot(4,1,4);plot(t,d-d1,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
title('(d)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
hold on;
plot(t1*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t2*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t3*[1,1],[-1,1],'r--','linewidth',1.5);
plot(t4*[1,1],[-1,1],'r--','linewidth',1.5);

annotation('textarrow',[0.6,0.66],[0.45,0.41],'color','r','linewidth',2,'Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','syn_dn2.eps');

yc_snr(d0,d)    %3.22
yc_snr(d0,d1)   %12.06

%% Frequency
d0_f=abs(fft(d0));d0_f=d0_f(1:128);
d_f=abs(fft(d));d_f=d_f(1:128);
d1_f=abs(fft(d1));d1_f=d1_f(1:128);
d2_f=abs(fft(d2));d2_f=d2_f(1:128);
d3_f=abs(fft(d3));d3_f=d3_f(1:128);
d4_f=abs(fft(d4));d4_f=d4_f(1:128);

dt=0.002;
Nf=128;df=1/dt/2/Nf;
f=[0:Nf-1]*df;
figure;plot(f,d0_f,'r','linewidth',2);


figure('units','normalized','Position',[0.2 0.4 0.4, 0.9]);
subplot(6,1,1);plot(f,d0_f,'r','linewidth',2);ylim([0,12]);xlim([0,250]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(6,1,2);plot(f,d_f,'r','linewidth',2);ylim([0,12]);xlim([0,250]);
ylabel('Amplitude','Fontsize',16);
title('(b)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);

 subplot(6,1,3);plot(f,d4_f,'r','linewidth',2);ylim([0,12]);xlim([0,250]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);

subplot(6,1,4);plot(f,d3_f,'r','linewidth',2);ylim([0,12]);xlim([0,250]);
ylabel('Amplitude','Fontsize',16);
title('(d)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
% annotation('textarrow',[0.6,0.66],[0.45,0.41],'color','r','linewidth',2,'Fontsize',20,'fontweight','bold');

 subplot(6,1,5);plot(f,d2_f,'r','linewidth',2);ylim([0,12]);xlim([0,250]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(e)','Fontsize',16);

subplot(6,1,6);plot(f,d1_f,'r','linewidth',2);ylim([0,12]);xlim([0,250]);
ylabel('Amplitude','Fontsize',16);
xlabel('Frequency (Hz)','Fontsize',16);
title('(f)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);

annotate;
print(gcf,'-depsc','-r300','syn_dn_f.eps');


%   figure('units','normalized','Position',[0.2 0.4 0.55, 0.35]);
figure('units','normalized');
imagesc(Gksvd0);colormap(jet);colorbar;caxis([-1,2.0]);
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_Gksvd0.eps');

figure('units','normalized');
imagesc(Gksvd);colormap(jet);caxis([-1,2.0]);colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_Gksvd.eps');

figure('units','normalized');
imagesc(D);colormap(jet);caxis([-0.5,0.5]);%colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Atom NO','Fontsize',16);
title('Dictionary Atoms','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_Dct.eps');

figure('units','normalized');
imagesc(Dksvd);colormap(jet);caxis([-0.5,0.5]);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Atom NO','Fontsize',16);
title('Dictionary Atoms','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_Dksvd.eps');

figure('units','normalized');
imagesc(X);colormap(jet);%caxis([-0.5,0.5]);%colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_X.eps');

figure('units','normalized');
imagesc(X1);colormap(jet);%caxis([-0.5,0.5]);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_X1.eps');


%% plot the first 64 atoms
figure('units','normalized','Position',[0.2 0.4 0.6, 0.6]);
for ia=1:16
    subplot(4,4,ia);plot([1:l1],dct(:,ia),'b','linewidth',2);
    set(gca,'Linewidth',1.5,'Fontsize',16);
         ytickformat('%.1f'); 
    if ismember(ia,[1,5,9,13])
      ylabel('Amplitude','Fontsize',16);  

                        else
        set(gca,'yticklabel',[]);
    end
   
        
    if ismember(ia,[13,14,15,16])
      xlabel('Sample NO','Fontsize',16);
                        else
        set(gca,'xticklabel',[]);
    end
    xlim([1,l1]);
end
print(gcf,'-depsc','-r300','syn_atoms.eps');

% ytickformat('%.1f');
figure('units','normalized','Position',[0.2 0.4 0.6, 0.6]);
for ia=1:16
    subplot(4,4,ia);plot([1:l1],Dksvd(:,ia),'b','linewidth',2);
    
    if ismember(ia,[4,5,6,7])
       subplot(4,4,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2); 
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    ytickformat('%.1f');
    if ismember(ia,[1,5,9,13])
      ylabel('Amplitude','Fontsize',16);  
      
                  else
        set(gca,'yticklabel',[]);
    end
   
        
    if ismember(ia,[13,14,15,16])
      xlabel('Sample NO','Fontsize',16);  
                  else
        set(gca,'xticklabel',[]);
    end
    xlim([1,l1]);
end
print(gcf,'-depsc','-r300','syn_atoms1.eps');



% test percentage
%case 1
% Gksvd=pthresh(Gksvd0,'ph',0.4);
XX1=Dksvd*Gksvd0;
dd1=yc_patch_inv(XX1,1,n1,n2,16,1,8,1);

XX2=Dksvd*yc_pthresh(Gksvd0,'ph',0.8);
dd2=yc_patch_inv(XX2,1,n1,n2,16,1,8,1);

XX3=Dksvd*yc_pthresh(Gksvd0,'ph',1.0);
dd3=yc_patch_inv(XX3,1,n1,n2,16,1,8,1);

XX4=Dksvd*yc_pthresh(Gksvd0,'ph',2.0);
dd4=yc_patch_inv(XX4,1,n1,n2,16,1,8,1);

 figure('units','normalized','Position',[0.2 0.4 0.4, 0.6]);
subplot(4,1,1);plot(t,dd2,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(4,1,2);plot(t,dd3,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);

subplot(4,1,3);plot(t,dd4,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);

subplot(4,1,4);plot(t,dd1,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
title('(d)','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
% annotation('textarrow',[0.6,0.66],[0.45,0.41],'color','r','linewidth',2,'Fontsize',20,'fontweight','bold');

print(gcf,'-depsc','-r300','syn_ddn.eps');



