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


% x=0:.1:2*pi;
% y=sin(x);
% figure;plot(x,y);
% 


clc;clear;close all;

%%
ne=20;%number of events;
dt = 4./1000;
tmax = dt*255;
h = [0:10:1280-10];
tau=linspace(0.1,1.8,ne);
v0=linspace(1800,3000,ne);%exact
randn('state',201819202122);amp = randn(ne,1);%[1., -1.,1];
f0 = 20;
snr = 200;%default snr=2
L = 6;
seed=201517;
dh=hevents(dt,f0,tmax,h,tau,v0,amp,snr,L,seed);dh=yc_scale(dh,2);
figure;imagesc(dh);

tau1=[0,0];
p=[0.004,-0.004];
amp1=[1,1];
dl = levents(dt,f0,tmax,h,tau1,p,amp1,snr,L,seed);dl(:,32:end)=0;dl=yc_scale(dl,2);
dl2 = levents(dt,f0,tmax,h,tau1,[0.002,-0.002],amp1,snr,L,seed);dl2=yc_scale(dl2,2);
% dl2=0;
randn('state',202021);
n=0.1*randn(size(dl));
% n=0;
d=dh+dl+dl2+n;
figure;imagesc(d);

figure('units','normalized','Position',[0.2 0.4 0.53, 1],'color','w');
imagesc(d);colormap(jet);caxis([-0.50,0.50]);%colormap(seis);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Trace','Fontsize',20,'fontweight','bold');
set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
% title('Synthetic Example','Fontsize',20,'fontweight','bold');
print(gcf,'-depsc','-r300','sta.eps');


% ne=20;%number of events;
% dt = 4./1000;
% tmax = dt*511;
% h = [-640:10:630];
% tau=linspace(0.1,1.8,ne);
% v0=linspace(1800,3000,ne);%exact
% randn('state',201819202122);amp = randn(ne,1);%[1., -1.,1];
% f0 = 20;
% snr = 2000;%default snr=2
% L = 6;
% seed=201517;
% dh=hevents(dt,f0,tmax,h,tau,v0,amp,snr,L,seed);dh=yc_scale(dh,2);
% 
% tau1=[0,0];
% p=[0.004,-0.004];
% amp1=[1,1];
% dl = levents(dt,f0,tmax,h,tau1,p,amp1,snr,L,seed);dl=yc_scale(dl,2);
% % dl = dl+levents(dt,f0,tmax,h,tau1,p*1.5,amp1,snr,L,seed);%+levents(dt,f0,tmax,h,tau1,p*2.5,amp1,snr,L,seed);;
% randn('state',202021);
% % n=0.1*randn(size(dl));
% n=0;
% d=dh+dl+n;



% [n1,n2]=size(d);
% t=[0:n1-1]*dt;
% figure;imagesc(h,t,d);colormap(seis);


d0=d;%clean

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
%[Dksvd,Gksvd,err] = ksvd(params,'');
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

%% atoms

% option 2
inds=[9,10,11,12,21,22,23,24,25,26,29,38,44,52,53,55,56,57,62];

inds=[9,10,18,21,22,23,24,25,26,29,38,44,52,53,55,56,57,62];%in backup statistics
% 
% inds=[9,10,15,16,25,26,29,30,38,46,53,55,57,58,60];

inds=[8,9,10,11,12,13,14,15,16,17,19,23,24,25,26,27,28,29,30,31,32,40,41,42,43,44,45,46,47,48,56,57,58,59,60,61,62,63,64];


figure('units','normalized','Position',[0.2 0.4 0.53, 1]);
for ia=1:64
    subplot(8,8,ia);imagesc(linspace(-1,1,c1),linspace(-1,1,c1),reshape(Dksvd(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
%     title(num2str(var(Dksvd(:,ia))));
    if ismember(ia,inds)
        draw_circle(0,0,0.5,'k');
    end
    set(gca,'Linewidth',1.5,'Fontsize',16);
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
end
print(gcf,'-depsc','-r300','sta_atoms.eps');

%% select atoms
%option 1
% inds=[11,12,13,14,15,19,20,21,22,23,25,28,29,38,39,41,43,44,46,49,50,52,55,56,60,61,63];
% Dksvd_o1=Dksvd;
% Dksvd_o1(:,inds)=0;
% X2_o1=Dksvd_o1*Gksvd_omp;
% d2_o1=yc_patch_inv(X2_o1,1,n1,n2,l1,l2,l1/2,l2/2);
% figure;imagesc([d,d2_o1,d-d2_o1]);colormap(seis);
% yc_snr(dh,d)
% yc_snr(dh,d2_o1)



m=abs(yc_mean(Dksvd));m=yc_scale(m);
md=abs(yc_median(Dksvd));md=yc_scale(md);
v=yc_var(Dksvd);v=yc_scale(v);
k=yc_kurtosis(Dksvd);k=yc_scale(k);
s=abs(yc_skewness(Dksvd));s=yc_scale(s);

figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
plot(1:64,m,'-mp','linewidth',2); hold on;
plot(1:64,md,'-gd','linewidth',2);
plot(1:64,v,'-c^','linewidth',2);
plot(1:64,s,'-bo','linewidth',2);
plot(1:64,k,'-rs','linewidth',2);

%% create windows
wins=inds(:);
wins=[wins-0.5,wins+0.5];
nwin=size(wins,1);
amp=1;
dt=1;
for iw=1:nwin

        ws=wins(iw,1);we=wins(iw,2);
        x=[ws,we,we,ws]*dt;y=[0,0,amp,amp];
        patch(x,y,'b','FaceAlpha',0.1);

end
legend('Mean','Median','Variance','Skewness','Kurtosis','location','best');
xlim([0,65]);
ylabel('Value','Fontsize',20,'fontweight','bold');
xlabel('Atom NO','Fontsize',20,'fontweight','bold');
% title('Noise','Fontsize',20,'fontweight','bold');
  set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
  print(gcf,'-depsc','-r300','sta_curves.eps');


b_man=zeros(64,1);
b_man(inds)=1;

% c_m=sum(m(:).*b_man);
% c_md=sum(md(:).*b_man);
% c_v=sum(v(:).*b_man);
% c_k=sum(k(:).*b_man);
% c_s=sum(s(:).*b_man);
% e_m=yc_rms(m(:)-b_man)
% e_md=yc_rms(md(:)-b_man)
% e_v=yc_rms(v(:)-b_man)
% e_s=yc_rms(s(:)-b_man)
% e_k=yc_rms(k(:)-b_man)
% 
% 
% [c_m,c_md,c_v,c_k,c_s]
% 
N=length(inds);
inds_m=sort(m,'descend');

[~,ii]=sort(m,'descend');
inds_m=ii(1:N);

[~,ii]=sort(md,'descend');
inds_md=ii(1:N);

[~,ii]=sort(v,'descend');
inds_v=ii(1:N);

[~,ii]=sort(s,'descend');
inds_s=ii(1:N);

[~,ii]=sort(k,'descend');
inds_k=ii(1:N);


b_m=zeros(64,1);
b_m(inds_m)=1;

b_md=zeros(64,1);
b_md(inds_md)=1;

b_v=zeros(64,1);
b_v(inds_v)=1;

b_s=zeros(64,1);
b_s(inds_s)=1;

b_k=zeros(64,1);
b_k(inds_k)=1;

figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
plot(1:64,b_m,'-mp','linewidth',2); hold on;
plot(1:64,b_md,'-gd','linewidth',2);
plot(1:64,b_v,'-c^','linewidth',2);
plot(1:64,b_s,'-bo','linewidth',2);
plot(1:64,b_k,'-rs','linewidth',2);
plot(1:64,b_man,'-k*','linewidth',2);
legend('Mean','Median','Variance','Skewness','Kurtosis','Manual','location','best');
xlim([0,65])
legend('Mean','Median','Variance','Skewness','Kurtosis','Solution','location','best');
ylabel('Result','Fontsize',20,'fontweight','bold');
xlabel('Atom NO','Fontsize',20,'fontweight','bold');
% title('Noise','Fontsize',20,'fontweight','bold');
  set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
  print(gcf,'-depsc','-r300','sta_result.eps');

  
sum(abs(b_man-b_m))
sum(abs(b_man-b_md))
sum(abs(b_man-b_v))
sum(abs(b_man-b_s))
sum(abs(b_man-b_k))


figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
plot(1:64,abs(b_man-b_m),'-mp','linewidth',2); hold on;
plot(1:64,abs(b_man-b_md),'-gd','linewidth',2);
plot(1:64,abs(b_man-b_v),'-c^','linewidth',2);
plot(1:64,abs(b_man-b_s),'-bo','linewidth',2);
plot(1:64,abs(b_man-b_k),'-rs','linewidth',2);
xlim([0,65])
legend('Mean (E=42)','Median (E=32)','Variance (E=22)','Skewness (E=18)','Kurtosis (E=6)','location','best');
ylabel('Error','Fontsize',20,'fontweight','bold');
xlabel('Atom NO','Fontsize',20,'fontweight','bold');
% title('Noise','Fontsize',20,'fontweight','bold');
  set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
  print(gcf,'-depsc','-r300','sta_err.eps');

% %% Noise
% 
% ne=20;%number of events;
% dt = 4./1000;
% tmax = dt*255;
% h = [0:10:1280-10];
% tau=linspace(0.1,1.8,ne);
% v0=linspace(1800,3000,ne);%exact
% randn('state',201819202122);amp = randn(ne,1);%[1., -1.,1];
% f0 = 20;
% snr = 200;%default snr=2
% L = 6;
% seed=201517;
% dh=hevents(dt,f0,tmax,h,tau,v0,amp,snr,L,seed);dh=yc_scale(dh,2);
% figure;imagesc(dh);
% 
% tau1=[0,0];
% p=[0.004,-0.004];
% amp1=[1,1];
% dl = levents(dt,f0,tmax,h,tau1,p,amp1,snr,L,seed);dl(:,32:end)=0;dl=yc_scale(dl,2);
% dl2 = levents(dt,f0,tmax,h,tau1,[0.002,-0.002],amp1*0.5,snr,L,seed);dl2=yc_scale(dl2,2);
% 
% randn('state',202021);
% n=0.1*randn(size(dl));
% 
% d=dh+dl+dl2+n;
% figure;imagesc(d);
% 
% % ne=20;%number of events;
% % dt = 4./1000;
% % tmax = dt*511;
% % h = [-640:10:630];
% % tau=linspace(0.1,1.8,ne);
% % v0=linspace(1800,3000,ne);%exact
% % randn('state',201819202122);amp = randn(ne,1);%[1., -1.,1];
% % f0 = 20;
% % snr = 2000;%default snr=2
% % L = 6;
% % seed=201517;
% % dh=hevents(dt,f0,tmax,h,tau,v0,amp,snr,L,seed);dh=yc_scale(dh,2);
% % 
% % tau1=[0,0];
% % p=[0.004,-0.004];
% % amp1=[1,1];
% % dl = levents(dt,f0,tmax,h,tau1,p,amp1,snr,L,seed);dl=yc_scale(dl,2);
% % % dl = dl+levents(dt,f0,tmax,h,tau1,p*1.5,amp1,snr,L,seed);%+levents(dt,f0,tmax,h,tau1,p*2.5,amp1,snr,L,seed);;
% % randn('state',202021);
% % % n=0.1*randn(size(dl));
% % n=0;
% % d=dh+dl+n;
% 
% 
% 
% % [n1,n2]=size(d);
% % t=[0:n1-1]*dt;
% % figure;imagesc(h,t,d);colormap(seis);
% 
% 
% d0=d;%clean
% 
% % patch size l1*l2
% l1=8;l2=8;
% 
% c1=8;c2=16;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
% %% DCT dictionary (dctmtx will generates orthogonal transform)
% dct=zeros(c1,c2);
% for k=0:1:c2-1
%     V=cos([0:1:c1-1]'*k*pi/c2);
%     if k>0
%         V=V-mean(V);
%     end
%     dct(:,k+1)=V/norm(V);
% end
% DCT=kron(dct,dct);%2D DCT dictionary (64,256)
% % 
% %% plot the first 64 atoms
% figure;
% for ia=1:64
%         subplot(8,8,ia);imagesc(reshape(DCT(:,ia),c1,c1));
% end
% 
% %% testing patch/patch_inv function
% % decompose the image into patches:
% X=yc_patch(d,1,l1,l2,l1/2,l2/2);
% % OMP for G
% nd=size(X,2);               %in this case = 3969
% G=zeros(c2*c2,1);
% K=3;
% tic
% for i2=1:nd
%     G(:,i2)=yc_omp0(DCT,X(:,i2),K);
%     if mod(i2,50)==0
%         fprintf('i2=%d/%d is finished\n',i2,nd);
%     end
% end
% toc
% X1=DCT*G;
% % 
% % insert patches into the image
% [n1,n2]=size(d);
% d1=yc_patch_inv(X1,1,n1,n2,l1,l2,l1/2,l2/2);
% figure;imagesc([d,d1,d-d1]);colormap(seis);
% yc_snr(d0,d1)
% % 
% % 
% %% K-SVD
% addpath(genpath('~/chenyk.data2/various/packages/toolbox/Dictionary'));
% 
% params.data = X;
% params.Tdata = 3;%sparsity level
% params.initdict=DCT;
% params.dictsize = c2.*c2;
% params.dictsize = 64;
% % params.dictsize = 16;
% params.codemode = 'sparsity';
% % params.codemode = 'error';
% params.iternum = 30;
% params.memusage = 'high';
% 
% seed=201819;
% rng(seed,'twister');
% randn('state',seed);
% rand('state',seed);
% [Dksvd,Gksvd,err] = ksvd(params,'');
% 
% % Gksvd0=Gksvd;
% % Gksvd=yc_pthresh(Gksvd0,'ph',1);
% X2=Dksvd*Gksvd;
% 
% for i2=1:nd
%     Gksvd_omp(:,i2)=yc_omp0(Dksvd,X(:,i2),K);
%     if mod(i2,50)==0
%         fprintf('i2=%d/%d is finished\n',i2,nd);
%     end
% end
% X2_omp=Dksvd*Gksvd_omp;
% 
% [n1,n2]=size(d);
% d2=yc_patch_inv(X2,1,n1,n2,l1,l2,l1/2,l2/2);
% 
% d2_omp=yc_patch_inv(X2_omp,1,n1,n2,l1,l2,l1/2,l2/2);
% 
% 
% figure;imagesc([d,d2,d-d2]);colormap(seis);
% figure;imagesc([d,d2_omp,d-d2_omp]);colormap(seis);
% yc_snr(d0,d2)
% 
% %% atoms
% figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
% for ia=1:64
%     subplot(8,8,ia);imagesc(reshape(Dksvd(:,ia),c1,c1));colormap(jet);caxis([-0.1,0.1]);
% %     title(num2str(var(Dksvd(:,ia))));
%     set(gca,'Linewidth',1.5,'Fontsize',16);
%     set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
% end
% 
% %% select atoms
% %option 1
% % inds=[11,12,13,14,15,19,20,21,22,23,25,28,29,38,39,41,43,44,46,49,50,52,55,56,60,61,63];
% % Dksvd_o1=Dksvd;
% % Dksvd_o1(:,inds)=0;
% % X2_o1=Dksvd_o1*Gksvd_omp;
% % d2_o1=yc_patch_inv(X2_o1,1,n1,n2,l1,l2,l1/2,l2/2);
% % figure;imagesc([d,d2_o1,d-d2_o1]);colormap(seis);
% % yc_snr(dh,d)
% % yc_snr(dh,d2_o1)
% 
% % option 2
% inds=[9,10,11,12,21,22,23,24,25,26,29,38,44,52,53,55,56,57,62];
% 
% inds=[9,10,18,21,22,23,24,25,26,29,38,44,52,53,55,56,57,62];
% 
% 
% m=abs(yc_mean(Dksvd));m=yc_scale(m);
% md=abs(yc_median(Dksvd));md=yc_scale(md);
% v=yc_var(Dksvd);v=yc_scale(v);
% k=yc_kurtosis(Dksvd);k=yc_scale(k);
% s=abs(yc_skewness(Dksvd));s=yc_scale(s);
% 
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% plot(1:64,m,'-mp','linewidth',2); hold on;
% plot(1:64,md,'-gd','linewidth',2);
% plot(1:64,v,'-c^','linewidth',2);
% plot(1:64,s,'-bo','linewidth',2);
% plot(1:64,k,'-rs','linewidth',2);
% 
% %% create windows
% wins=inds(:);
% wins=[wins-0.5,wins+0.5];
% nwin=size(wins,1);
% amp=1;
% dt=1;
% for iw=1:nwin
% 
%         ws=wins(iw,1);we=wins(iw,2);
%         x=[ws,we,we,ws]*dt;y=[0,0,amp,amp];
%         patch(x,y,'b','FaceAlpha',0.1);
% 
% end
% legend('Mean','Median','Variance','Skewness','Kurtosis','location','best');
% xlim([0,65]);
% ylabel('Value','Fontsize',20,'fontweight','bold');
% xlabel('Atom N','Fontsize',20,'fontweight','bold');
% % title('Noise','Fontsize',20,'fontweight','bold');
%   set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
%   print(gcf,'-depsc','-r300','sta_result_n.eps');
% 
% 
% b_man=zeros(64,1);
% b_man(inds)=1;
% 
% % c_m=sum(m(:).*b_man);
% % c_md=sum(md(:).*b_man);
% % c_v=sum(v(:).*b_man);
% % c_k=sum(k(:).*b_man);
% % c_s=sum(s(:).*b_man);
% % e_m=yc_rms(m(:)-b_man)
% % e_md=yc_rms(md(:)-b_man)
% % e_v=yc_rms(v(:)-b_man)
% % e_s=yc_rms(s(:)-b_man)
% % e_k=yc_rms(k(:)-b_man)
% % 
% % 
% % [c_m,c_md,c_v,c_k,c_s]
% % 
% N=length(inds);
% inds_m=sort(m,'descend');
% 
% [~,ii]=sort(m,'descend');
% inds_m=ii(1:N);
% 
% [~,ii]=sort(md,'descend');
% inds_md=ii(1:N);
% 
% [~,ii]=sort(v,'descend');
% inds_v=ii(1:N);
% 
% [~,ii]=sort(s,'descend');
% inds_s=ii(1:N);
% 
% [~,ii]=sort(k,'descend');
% inds_k=ii(1:N);
% 
% 
% b_m=zeros(64,1);
% b_m(inds_m)=1;
% 
% b_md=zeros(64,1);
% b_md(inds_md)=1;
% 
% b_v=zeros(64,1);
% b_v(inds_v)=1;
% 
% b_s=zeros(64,1);
% b_s(inds_s)=1;
% 
% b_k=zeros(64,1);
% b_k(inds_k)=1;
% 
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% plot(1:64,b_m,'-mp','linewidth',2); hold on;
% plot(1:64,b_md,'-gd','linewidth',2);
% plot(1:64,b_v,'-c^','linewidth',2);
% plot(1:64,b_s,'-bo','linewidth',2);
% plot(1:64,b_k,'-rs','linewidth',2);
% plot(1:64,b_man,'-k*','linewidth',2);
% legend('Mean','Median','Variance','Skewness','Kurtosis','Manual','location','best');
% xlim([0,65])
% legend('Mean','Median','Variance','Skewness','Kurtosis','location','best');
% ylabel('Result','Fontsize',20,'fontweight','bold');
% xlabel('Atom N','Fontsize',20,'fontweight','bold');
% % title('Noise','Fontsize',20,'fontweight','bold');
%   set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
%   print(gcf,'-depsc','-r300','sta_result_n.eps');
% 
%   
% sum(abs(b_man-b_m))
% sum(abs(b_man-b_md))
% sum(abs(b_man-b_v))
% sum(abs(b_man-b_s))
% sum(abs(b_man-b_k))
% 
% 
% figure('units','normalized','Position',[0.2 0.4 0.8, 0.6],'color','w');
% plot(1:64,abs(b_man-b_m),'-mp','linewidth',2); hold on;
% plot(1:64,abs(b_man-b_md),'-gd','linewidth',2);
% plot(1:64,abs(b_man-b_v),'-c^','linewidth',2);
% plot(1:64,abs(b_man-b_s),'-bo','linewidth',2);
% plot(1:64,abs(b_man-b_k),'-rs','linewidth',2);
% xlim([0,65])
% legend('Mean (E=42)','Median (E=32)','Variance (E=22)','Skewness (E=18)','Kurtosis (E=6)','location','best');
% ylabel('Error','Fontsize',20,'fontweight','bold');
% xlabel('Atom N','Fontsize',20,'fontweight','bold');
% % title('Noise','Fontsize',20,'fontweight','bold');
%   set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
%   print(gcf,'-depsc','-r300','sta_err_n.eps');
% 
  
 