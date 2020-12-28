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

% DATA Source:
% http://certmapper.cr.usgs.gov/data/NPRA/SEISMIC/1981/31_81/PROCESSED/31_81_PR.SGY

fid=fopen('alaska_1501_534_4ms.bin','r');
d=fread(fid,[1501,534],'float');
d=yc_scale(d,2);
figure;imagesc(d);caxis([-0.2,0.2]);colormap(seis);

d2=fx_mssa(d,1,100,0.004,10,0);
figure;imagesc([d,d2,d-d2]);caxis([-0.1,0.1]);colormap(gray);

%%

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
G=zeros(c2*c2,1);%d1=d2_o1;

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

seed=201819;%fine
seed=2018192021;%some energy
seed=0123456789;
rng(seed,'twister');
randn('state',seed);
rand('state',seed);
%[Dksvd,Gksvd,err] = ksvd(params,'');
[Dksvd,Gksvd] = yc_ksvd(X,param);
norm(Dksvd)

% Gksvd0=Gksvd;
% Gksvd=yc_pthresh(Gksvd0,'ph',1);
X2=Dksvd*Gksvd;

% for i2=1:nd
%     Gksvd_omp(:,i2)=yc_omp0(Dksvd,X(:,i2),K);
%     if mod(i2,50)==0
%         fprintf('i2=%d/%d is finished\n',i2,nd);
%     end
% end
% X2_omp=Dksvd*Gksvd_omp;
% d2_omp=yc_patch_inv(X2_omp,1,n1,n2,l1,l2,l1/2,l2/2);
% figure;imagesc([d,d2_omp,d-d2_omp]);colormap(seis);

[n1,n2]=size(d);
d2=yc_patch_inv(X2,1,n1,n2,l1,l2,l1/2,l2/2);



figure;imagesc([d,d2,d-d2]);colormap(seis);


%% atoms
figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
for ia=1:64
    subplot(8,8,ia);imagesc(reshape(Dksvd(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end

%% coherent noise removal: select atoms (test different statistics)
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
% perc=35;%works fine
% perc=50;%works perfectly? I think so. %in original script
perc=55;
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
Dksvd_o1=Dksvd;
Dksvd_o1(:,inds)=0;

figure('units','normalized','Position',[0.2 0.4 0.4, 1.2]);
for ia=1:64
    subplot(16,8,ia);imagesc(reshape(Dksvd(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
for ia=1:64
    subplot(16,8,ia+64);imagesc(reshape(Dksvd_o1(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
annotation(gcf,'line',[0.003472 0.994791],...
    [0.5169 0.5169],'Color',[1 0 0],'LineWidth',3);


X2_o1=Dksvd_o1*Gksvd;
d2_o1=yc_patch_inv(X2_o1,1,n1,n2,l1,l2,l1/2,l2/2);
d1=d2_o1;
% figure;imagesc([d,d2_o1,d-d2_o1]);caxis([-0.25,0.25]);colormap(seis);
figure;imagesc([d,d2_o1,d-d2_o1]);caxis([-0.10,0.10]);colormap(gray);

%Annotate
annotation(gcf,'rectangle',...
    [0.754016323633783 0.631578947368421 0.150880766501065 0.0973684210526314],...
    'Color',[1 0 0],...
    'LineWidth',3);
annotation(gcf,'rectangle',...
    [0.732963692054835 0.249516441005803 0.0657205184714807 0.148676575384302],...
    'Color',[1 0 0],...
    'LineWidth',3);
annotation(gcf,'rectangle',...
    [0.868421052631579 0.129593810444874 0.0381578947368421 0.305609284332689],...
    'Color',[1 0 0],...
    'LineWidth',3);
annotation(gcf,'rectangle',...
    [0.648684210526316 0.439071566731141 0.0921052631578949 0.0947775628626705],...
    'Color',[1 0 0],...
    'LineWidth',3);
% Note that to remove coherent noise, the signal damages are inevitable


%% plot result
[n1,n2]=size(d);
dt=0.004;
t=[0:n1-1]*dt;
x=[1:n2];

  figure('units','normalized','Position',[0.2 0.4 1.0, 0.8],'color','w');
subplot(1,3,1);
imagesc(x,t,d);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) Alaska field data','Fontsize',20,'fontweight','bold');

subplot(1,3,2);
imagesc(x,t,d2_o1);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) Denoised data','Fontsize',20,'fontweight','bold');

subplot(1,3,3);
imagesc(x,t,d-d2_o1);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) Removed noise','Fontsize',20,'fontweight','bold');
real_anno;
% annotation(gcf,'rectangle',...
%     [0.761111111111112 0.232044198895028 0.0576388888888882 0.151012891344384],...
%     'Color',[1 0 0],...
%     'LineWidth',4);
% 
% % 创建 rectangle
% annotation(gcf,'rectangle',...
%     [0.8375 0.61878453038674 0.0659722222222222 0.0810313075506464],...
%     'Color',[1 0 0],...
%     'LineWidth',4);
% 
% % 创建 rectangle
% annotation(gcf,'rectangle',...
%     [0.871527777777778 0.186003683241252 0.0319444444444448 0.276243093922653],...
%     'Color',[1 0 0],...
%     'LineWidth',4);
% 
% % 创建 rectangle
% annotation(gcf,'rectangle',...
%     [0.693750000000001 0.421731123388582 0.068055555555555 0.114180478821363],...
%     'Color',[1 0 0],...
%     'LineWidth',4);
print(gcf,'-depsc','-r300','real.eps');

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
%% random noise removal
% Gksvd0=Gksvd;
% Gksvd=yc_pthresh(Gksvd0,'ph',4);
X3=Dksvd*Gksvd;
d3=yc_patch_inv(X3,1,n1,n2,l1,l2,l1/2,l2/2);
figure;imagesc([d,d3,d-d3]);caxis([-0.10,0.10]);colormap(gray);

%% fx EMD
d4 = fxemd(d,1,100,0.004,1,0);
figure;imagesc([d,d4,d-d4]);caxis([-0.10,0.10]);colormap(gray);

%% FK
d5= yc_fkt(d,'ps',40);
figure;imagesc([d,d5,d-d5]);caxis([-0.10,0.10]);colormap(gray);

%% comp
figure;imagesc([d2_o1,d2,d3,d4,d5]);caxis([-0.10,0.10]);colormap(gray);



figure('units','normalized','Position',[0.2 0.4 1.0, 0.8],'color','w');
subplot(1,3,3);
imagesc(x,t,d3);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(1,3,1);
imagesc(x,t,d2);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(1,3,2);
imagesc(x,t,d2);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');
real_anno;
print(gcf,'-depsc','-r300','real_dn.eps');


figure('units','normalized','Position',[0.2 0.4 1.0, 0.8],'color','w');
subplot(1,3,3);
imagesc(x,t,d-d3);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(1,3,1);
imagesc(x,t,d-d2);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(1,3,2);
imagesc(x,t,d-d4);colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');
real_anno;
print(gcf,'-depsc','-r300','real_n.eps');


%% zoom in 1
ind1=find(t>1.5 & t<2);
ind2=find(x>400);
figure('units','normalized','Position',[0.2 0.4 0.7, 1.0],'color','w');
subplot(2,3,1);
imagesc(x(ind2),t(ind1),d(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) Field data','Fontsize',20,'fontweight','bold');

subplot(2,3,2);
imagesc(x(ind2),t(ind1),d2(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) SSA','Fontsize',20,'fontweight','bold');

subplot(2,3,3);
imagesc(x(ind2),t(ind1),d4(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,3,5);
imagesc(x(ind2),t(ind1),d3(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) DL','Fontsize',20,'fontweight','bold');

subplot(2,3,6);
imagesc(x(ind2),t(ind1),d1(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(e) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_dn_z.eps');


nn1=d-d1;
nn2=d-d2;
n3=d-d3;
n4=d-d4;

figure('units','normalized','Position',[0.2 0.4 0.45, 1.0],'color','w');
subplot(2,2,1);
imagesc(x(ind2),t(ind1),nn2(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(2,2,2);
imagesc(x(ind2),t(ind1),n4(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,2,3);
imagesc(x(ind2),t(ind1),n3(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(2,2,4);
imagesc(x(ind2),t(ind1),nn1(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_n_z.eps');


%% Compact comparison
% figure('units','normalized','Position',[0.2 0.4 1.0, 1.0],'color','w');
% subplot(2,4,1);
% imagesc(x(ind2),t(ind1),d2(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(a) SSA','Fontsize',20,'fontweight','bold');
% 
% subplot(2,4,2);
% imagesc(x(ind2),t(ind1),d4(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% % set(gca,'yticklabel',[]);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(b) FX EMD','Fontsize',20,'fontweight','bold');
% 
% subplot(2,4,3);
% imagesc(x(ind2),t(ind1),d3(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(c) DL','Fontsize',20,'fontweight','bold');
% 
% subplot(2,4,4);
% imagesc(x(ind2),t(ind1),d1(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% % set(gca,'yticklabel',[]);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(d) Proposed','Fontsize',20,'fontweight','bold');
% 
% subplot(2,4,5);
% imagesc(x(ind2),t(ind1),nn2(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(a) SSA','Fontsize',20,'fontweight','bold');
% 
% subplot(2,4,6);
% imagesc(x(ind2),t(ind1),n4(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% % set(gca,'yticklabel',[]);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(b) FX EMD','Fontsize',20,'fontweight','bold');
% 
% subplot(2,4,7);
% imagesc(x(ind2),t(ind1),n3(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(c) DL','Fontsize',20,'fontweight','bold');
% 
% subplot(2,4,8);
% imagesc(x(ind2),t(ind1),nn1(ind1,ind2));colormap(seis);caxis([-0.15,0.15]);colormap(gray);
% % set(gca,'yticklabel',[]);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('(d) Proposed','Fontsize',20,'fontweight','bold');
% print(gcf,'-depsc','-r300','real_comp_z.eps');


%% zoom in 2
ind1=find(t>4 & t<5);
ind2=find(x>200 & x<300);
figure('units','normalized','Position',[0.2 0.4 0.7, 1.0],'color','w');
subplot(2,3,1);
imagesc(x(ind2),t(ind1),d(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) Field data','Fontsize',20,'fontweight','bold');

subplot(2,3,2);
imagesc(x(ind2),t(ind1),d2(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) SSA','Fontsize',20,'fontweight','bold');

subplot(2,3,3);
imagesc(x(ind2),t(ind1),d4(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,3,5);
imagesc(x(ind2),t(ind1),d3(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) DL','Fontsize',20,'fontweight','bold');

subplot(2,3,6);
imagesc(x(ind2),t(ind1),d1(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(e) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_dn_z2.eps');


nn1=d-d1;
nn2=d-d2;
n3=d-d3;
n4=d-d4;

figure('units','normalized','Position',[0.2 0.4 0.45, 1.0],'color','w');
subplot(2,2,1);
imagesc(x(ind2),t(ind1),nn2(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(2,2,2);
imagesc(x(ind2),t(ind1),n4(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,2,3);
imagesc(x(ind2),t(ind1),n3(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(2,2,4);
imagesc(x(ind2),t(ind1),nn1(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_n_z2.eps');


%% zoom in 3
ind1=find(t>3.4 & t<5.4);
ind2=find(x>450 );
figure('units','normalized','Position',[0.2 0.4 0.7, 1.0],'color','w');
subplot(2,3,1);
imagesc(x(ind2),t(ind1),d(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) Field data','Fontsize',20,'fontweight','bold');

subplot(2,3,2);
imagesc(x(ind2),t(ind1),d2(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) SSA','Fontsize',20,'fontweight','bold');

subplot(2,3,3);
imagesc(x(ind2),t(ind1),d4(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,3,5);
imagesc(x(ind2),t(ind1),d3(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) DL','Fontsize',20,'fontweight','bold');

subplot(2,3,6);
imagesc(x(ind2),t(ind1),d1(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(e) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_dn_z3.eps');


nn1=d-d1;
nn2=d-d2;
n3=d-d3;
n4=d-d4;

figure('units','normalized','Position',[0.2 0.4 0.45, 1.0],'color','w');
subplot(2,2,1);
imagesc(x(ind2),t(ind1),nn2(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(2,2,2);
imagesc(x(ind2),t(ind1),n4(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,2,3);
imagesc(x(ind2),t(ind1),n3(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(2,2,4);
imagesc(x(ind2),t(ind1),nn1(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_n_z3.eps');



%% zoom in 4
ind1=find(t>2.8 & t<3.8);
ind2=find(x>20 & x<150 );
figure('units','normalized','Position',[0.2 0.4 0.7, 1.0],'color','w');
subplot(2,3,1);
imagesc(x(ind2),t(ind1),d(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) Field data','Fontsize',20,'fontweight','bold');

subplot(2,3,2);
imagesc(x(ind2),t(ind1),d2(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) SSA','Fontsize',20,'fontweight','bold');

subplot(2,3,3);
imagesc(x(ind2),t(ind1),d4(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,3,5);
imagesc(x(ind2),t(ind1),d3(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) DL','Fontsize',20,'fontweight','bold');

subplot(2,3,6);
imagesc(x(ind2),t(ind1),d1(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(e) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_dn_z4.eps');


nn1=d-d1;
nn2=d-d2;
n3=d-d3;
n4=d-d4;

figure('units','normalized','Position',[0.2 0.4 0.45, 1.0],'color','w');
subplot(2,2,1);
imagesc(x(ind2),t(ind1),nn2(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(a) SSA','Fontsize',20,'fontweight','bold');

subplot(2,2,2);
imagesc(x(ind2),t(ind1),n4(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(b) FX EMD','Fontsize',20,'fontweight','bold');

subplot(2,2,3);
imagesc(x(ind2),t(ind1),n3(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(c) DL','Fontsize',20,'fontweight','bold');

subplot(2,2,4);
imagesc(x(ind2),t(ind1),nn1(ind1,ind2));colormap(seis);caxis([-0.10,0.10]);colormap(gray);
% set(gca,'yticklabel',[]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
title('(d) Proposed','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','real_n_z4.eps');

%% plot the atoms
k=kurtosis(Dksvd);
r=yc_scale(k-min(k)); %radius

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(reshape(DCT(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','real_atom1.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(reshape(Dksvd(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','real_atom2.eps');

figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(reshape(Dksvd_o1(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','real_atom3.eps');

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
print(gcf,'-depsc','-r300','real_atom4.eps');

%% for fun
figure('units','normalized','Position',[0.2 0.4 0.53, 1]);%carefully adjust so that each box is a square
for ia=1:64
    subplot(8,8,ia);
    if ismember(ia,inds)
    plot(0, 0, '.g', 'MarkerSize',r(ia)*190+.1);    
    else
    plot(0, 0, '.r', 'MarkerSize',r(ia)*190+.1);
    end
%     title(num2str(var(Dksvd(:,ia))));
    set(gca,'Linewidth',1.5,'Fontsize',16,'XColor',[204 204 204]/255,'YColor',[204 204 204]/255,'Color','k');
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    
    set(gca,'Linewidth',1.5,'Fontsize',16,'XColor','k','YColor','k','Color','k');
end
set(gcf,'Color','k')
axis off;

set(gcf, 'InvertHardCopy', 'off'); 
print(gcf,'-depsc','-r300','atoms_for_fun.eps');
% annotation(gcf,'line',[0.003472 0.994791],...
%     [0.5169 0.5169],'Color',[1 0 0],'LineWidth',3);
% 


