clc;clear;close all;
% Demo for Statistics-guided dictionary learning for automatic coherent noise suppression
% Prepared By Yatong Zhou, Jian Yang, Hang Wang, Guangtan and Yangkang Chen
% Dec, 27, 2020
% 
% NOTE that the initial code of Zhou et al. (2021) uses the standard KSVD toolbox
% However, this demo script uses the standalone yc_ksvd.m function.
%
% Key Reference
% Zhou, Y., J. Yang, H. Wang, G. Huang, and Y. Chen, 2021, Statistics-guided dictionary learning for automatic coherent noise suppression, IEEE Transactions on Geoscience and Remote Sensing, doi: 10.1109/TGRS.2020.3039738.
% 
% For more details about dictionary learning and denoising, please refer to
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51-KS61.
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

%%
ne=20;%number of events;
dt = 4./1000;
tmax = dt*511;
h = [-640:10:630];
tau=linspace(0.1,1.8,ne);
v0=linspace(1800,3000,ne);%exact
randn('state',201819202122);amp = randn(ne,1);%[1., -1.,1];
f0 = 20;
snr = 2000;%default snr=2
L = 6;
seed=201517;
dh=hevents(dt,f0,tmax,h,tau,v0,amp,snr,L,seed);dh=yc_scale(dh,2);

tau1=[0,0];
p=[0.004,-0.004];
amp1=[1,1];
dl = levents(dt,f0,tmax,h,tau1,p,amp1,snr,L,seed);dl=yc_scale(dl,2);
% dl = dl+levents(dt,f0,tmax,h,tau1,p*1.5,amp1,snr,L,seed);%+levents(dt,f0,tmax,h,tau1,p*2.5,amp1,snr,L,seed);;
randn('state',202021);
n=0.1*randn(size(dl));
d=dh+dl+n;

[n1,n2]=size(d);
t=[0:n1-1]*dt;x=1:n2;
figure;imagesc(h,t,d);colormap(seis);
d0=dh;%clean

figure('units','normalized','Position',[0.2 0.4 1.0, 0.8],'color','w');
subplot(1,3,1);
imagesc(x,t,dh);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) Hyperbolic signals','Fontsize',20,'fontweight','bold');

subplot(1,3,2);
imagesc(x,t,dl);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) Coherent noise','Fontsize',20,'fontweight','bold');

subplot(1,3,3);
imagesc(x,t,d);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) Synthetic data','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','syn.eps');


%% Denoising starts here
% patch size l1*l2
l1=8;l2=8;
c1=8;c2=16;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
%% DCT dictionary (dctmtx will generates orthogonal transform)
dct=zeros(c1,c2);
for k=0:1:c2-1
    V=cos([0:1:c1-1]'*k*pi/c2);
    if k>0
        V=V-mean(V);
    end
    dct(:,k+1)=V/norm(V);
end
DCT=kron(dct,dct);%2D DCT dictionary (64,256)
% 
%% plot the first 64 atoms
figure;
for ia=1:64
        subplot(8,8,ia);imagesc(reshape(DCT(:,ia),c1,c1));
end

%% testing patch/patch_inv function
% decompose the image into patches:
X=yc_patch(d,1,l1,l2,l1/2,l2/2);
% OMP for G
nd=size(X,2);               %in this case = 3969
G=zeros(c2*c2,1);
K=3;
tic
for i2=1:nd
    G(:,i2)=yc_omp0(DCT,X(:,i2),K);
    if mod(i2,50)==0
        fprintf('i2=%d/%d is finished\n',i2,nd);
    end
end
toc
X1=DCT*G;
% 
% insert patches into the image
[n1,n2]=size(d);
d1=yc_patch_inv(X1,1,n1,n2,l1,l2,l1/2,l2/2);
figure;imagesc([d,d1,d-d1]);colormap(seis);
yc_snr(d0,d1)
% 
% 
%% K-SVD
%NOTE that the initial code uses the standard KSVD toolbox
%However, this demo script uses the standalone yc_ksvd.m function.
param.T = 3;%sparsity level
param.D=DCT;
param.niter=30;
param.K=64;     %number of atoms, dictionary size
natom=param.K;
param.mode=1;   %1: sparsity; 0: error

seed=201819;
rng(seed,'twister');
randn('state',seed);
rand('state',seed);
% [Dksvd,Gksvd,err] = ksvd(params,'');
[Dksvd,Gksvd] = yc_ksvd(X,param);

% Gksvd0=Gksvd;
% Gksvd=yc_pthresh(Gksvd0,'ph',1);
X2=Dksvd*Gksvd;

for i2=1:nd
    Gksvd_omp(:,i2)=yc_omp0(Dksvd,X(:,i2),K);
    if mod(i2,50)==0
        fprintf('i2=%d/%d is finished\n',i2,nd);
    end
end
X2_omp=Dksvd*Gksvd_omp;

[n1,n2]=size(d);
d2=yc_patch_inv(X2,1,n1,n2,l1,l2,l1/2,l2/2);

d2_omp=yc_patch_inv(X2_omp,1,n1,n2,l1,l2,l1/2,l2/2);


figure;imagesc([d,d2,d-d2]);colormap(seis);
figure;imagesc([d,d2_omp,d-d2_omp]);colormap(seis);
yc_snr(d0,d2)



%% select atoms (test different statistics)
%option 1
% inds=[11,12,13,14,15,19,20,21,22,23,25,28,29,38,39,41,43,44,46,49,50,52,55,56,60,61,63]; %17.0057 (human eye is worse than the kurtosis)

% vars=var(Dksvd); %13.9410 N=21
% [~,ii]=sort(vars,'descend');
% inds=ii(1:21);
% 
k=kurtosis(Dksvd);%16.6513 N=21; 16.9889 N=22; 17.1545 N=23; 17.3897 N=24; 17.6918 N=25; 17.1196 N=26
[ks,ii]=sort(k,'descend');
% figure;stem(ks);
% ks_ratio=[ks(1:end-1)./ks(2:end)];% ks_dif=diff(ks);figure;stem(ks_dif);
% figure;stem(ks_ratio);

perc=35;%works fine
perc=50;%works fine
perc=75;%works perfectly
perc=74;%works perfectly
% perc=39;%works 
tt=round(natom*(100-perc)/100);
inds=ii(1:tt);
% 
% m=mean(Dksvd);%3.1871 N=21
% [~,ii]=sort(m,'descend');
% inds=ii(1:21);
% % 
% s=skewness(Dksvd);%11.1845 N=21
% [~,ii]=sort(s,'descend');
% inds=ii(1:21);
%
% comp
% mean,median,variance,skewness,kurtosis

Dksvd_o1=Dksvd;
Dksvd_o1(:,inds)=0;

k=kurtosis(Dksvd);
r=yc_scale(k-min(k)); %radius

%% atoms
figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(reshape(DCT(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','syn_atom1.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(reshape(Dksvd(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','syn_atom2.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(reshape(Dksvd_o1(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','syn_atom3.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);%carefully adjust so that each box is a square
for ia=1:64
    subplot(8,8,ia);
    if ismember(ia,inds)
    plot(0, 0, '.g', 'MarkerSize',r(ia)*190+.1);    
    else
    plot(0, 0, '.r', 'MarkerSize',r(ia)*190+.1);
    end
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','syn_atom4.eps');





X2_o1=Dksvd_o1*Gksvd_omp;
d2_o1=yc_patch_inv(X2_o1,1,n1,n2,l1,l2,l1/2,l2/2);

X2=Dksvd*Gksvd_omp;
d3=yc_patch_inv(X2,1,n1,n2,l1,l2,l1/2,l2/2);

figure;imagesc([d,d2_o1,d-d2_o1]);colormap(seis);
yc_snr(dh,d)
yc_snr(dh,d2_o1)

%% Patches
figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc(X);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Input Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
annotation(gcf,'rectangle',...
    [0.448916666666667 0.108256880733945 0.0378888888888889 0.81651376146789],...
    'Color',[1 0 0],...
    'LineWidth',4);
print(gcf,'-depsc','-r300','syn_X.eps');

figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc(X2_o1);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Output Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
annotation(gcf,'rectangle',...
    [0.448916666666667 0.108256880733945 0.0378888888888889 0.81651376146789],...
    'Color',[1 0 0],...
    'LineWidth',4);
print(gcf,'-depsc','-r300','syn_X1.eps');

figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc([X-X2_o1]);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Noise patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
annotation(gcf,'rectangle',...
    [0.448916666666667 0.108256880733945 0.0378888888888889 0.81651376146789],...
    'Color',[1 0 0],...
    'LineWidth',4);
print(gcf,'-depsc','-r300','syn_X1_n.eps');

figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc(X2);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Output Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_X2.eps');

figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc([X-X2]);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Noise Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_X2_n.eps');

patch_t=1:64;
patch_x=1:size(X,2);
patch_inds=1625:1800;
figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc(patch_x(patch_inds),patch_t,[X(:,patch_inds)]);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Input Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);

annotation(gcf,'textarrow',[0.288194444444444 0.257638888888889],...
    [0.318709276231012 0.231797188318925],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Coherent noise'},...
    'LineWidth',4,...
    'FontWeight','bold',...
    'FontSize',30);

% 创建 textarrow
annotation(gcf,'textarrow',[0.292361111111111 0.261805555555556],...
    [0.652701371092673 0.565789283180587],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Coherent noise'},...
    'LineWidth',4,...
    'FontWeight','bold',...
    'FontSize',30);

% 创建 textarrow
annotation(gcf,'textarrow',[0.341666666666667 0.311111111111111],...
    [0.840448406665794 0.753536318753707],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Coherent noise'},...
    'LineWidth',4,...
    'FontWeight','bold',...
    'FontSize',30);

print(gcf,'-depsc','-r300','syn_X_z.eps');

figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc(patch_x(patch_inds),patch_t,[X2_o1(:,patch_inds)]);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Output Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_X1_z.eps');

figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc(patch_x(patch_inds),patch_t,[X(:,patch_inds)-X2_o1(:,patch_inds)]);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Output Patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
annotation(gcf,'textarrow',[0.288194444444444 0.257638888888889],...
    [0.318709276231012 0.231797188318925],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Coherent noise'},...
    'LineWidth',4,...
    'FontWeight','bold',...
    'FontSize',30);

% 创建 textarrow
annotation(gcf,'textarrow',[0.292361111111111 0.261805555555556],...
    [0.652701371092673 0.565789283180587],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Coherent noise'},...
    'LineWidth',4,...
    'FontWeight','bold',...
    'FontSize',30);

% 创建 textarrow
annotation(gcf,'textarrow',[0.341666666666667 0.311111111111111],...
    [0.840448406665794 0.753536318753707],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Coherent noise'},...
    'LineWidth',4,...
    'FontWeight','bold',...
    'FontSize',30);
print(gcf,'-depsc','-r300','syn_X1_n_z.eps');

figure('units','normalized','Position',[0.2 0.4 0.6, 0.5],'color','w');
imagesc(patch_x(patch_inds),patch_t,[X(:,patch_inds)-X2(:,patch_inds)]);colormap(jet);caxis([-0.5,0.5]);colormap(seis);%caxis([0,2.0]);colorbar;
ylabel('Sample NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Noise patches','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
print(gcf,'-depsc','-r300','syn_X2_n_z.eps');


% [n1,n2]=size(d);
% dt=0.004;
% t=[0:n1-1]*dt;
% x=[1:n2];

d1=d2_o1;
d4 = fxemd(d,1,100,0.004,1,0);
figure;imagesc([d,d4,d-d4]);caxis([-0.50,0.50]);colormap(seis);

d2=fx_mssa(d,1,100,0.004,20,0);
figure;imagesc([d,d2,d-d2]);caxis([-0.50,0.50]);colormap(seis);



figure('units','normalized','Position',[0.2 0.4 1.0, 0.8],'color','w');
subplot(1,4,1);
imagesc(x,t,d2);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(1,4,2);
imagesc(x,t,d4);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');

subplot(1,4,3);
imagesc(x,t,d3);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(1,4,4);
imagesc(x,t,d1);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','syn_dn.eps');


figure('units','normalized','Position',[0.2 0.4 1.0, 0.8],'color','w');
subplot(1,4,1);
imagesc(x,t,d-d2);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(1,4,2);
imagesc(x,t,d-d4);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');

subplot(1,4,3);
imagesc(x,t,d-d3);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(1,4,4);
imagesc(x,t,d-d1);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','syn_n.eps');


%% SNRs
% Noisy: 3.86 dB
% SSA: 5.75 dB
% FXEMD: 8.18 dB
% DL: 9.55 dB
% Proposed: 13.61 dB


%% test my OMP and the state-of-the-art OMP
% figure;imagesc([Gksvd,Gksvd_omp,Gksvd-Gksvd_omp]);caxis([-0.1,0.1]);

%% localortho
% rect = zeros(3, 1);
% rect(1) = 10;
% rect(2) = 10;
% rect(3) = 1;
% %eps=0.00000000000001;%too strong, do not know why
% eps=0;
% niter=20;
% verb=1;
% [d2_ortho,n2_ortho,low]=yc_localortho(d2_o1,d-d2_o1,rect,niter,eps,verb);
% yc_snr(dh,d2_ortho)
% figure;imagesc([d,d2_ortho,d-d2_ortho]);colormap(seis);
% 
% figure;imagesc([d-d2_o1,n2_ortho]);colormap(seis);
% 

%% Different percentage

ps=[50,70,75,80];
nps=length(ps);
As={'(a)','(b)','(c)','(d)'};
Bs={'(e)','(f)','(g)','(h)'};
figure('units','normalized','Position',[0.2 0.4 0.7, 1.0],'color','w');
for ips=1:nps
tt_iter=round(natom*(100-ps(ips))/100);
inds_iter=ii(1:tt_iter);

Dksvd_o1_iter=Dksvd;
Dksvd_o1_iter(:,inds_iter)=0;

X2_iter=Dksvd_o1_iter*Gksvd_omp;
d1_iter=yc_patch_inv(X2_iter,1,n1,n2,l1,l2,l1/2,l2/2);

subplot(2,nps,ips);
imagesc(x,t,d1_iter);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title(strcat(As{ips},' p=',num2str(ps(ips)),'%'),'Fontsize',20,'fontweight','bold');

subplot(2,nps,ips+nps);
imagesc(x,t,d-d1_iter);colormap(seis);caxis([-0.50,0.50]);colormap(seis);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title(strcat(Bs{ips},' p=',num2str(ps(ips)),'%'),'Fontsize',20,'fontweight','bold');

fprintf('ps=%g,SNR=%g\n',ps(ips),yc_snr(dh,d1_iter));
end
print(gcf,'-depsc','-r300','syn_ps.eps');
