% Demo for microseismic data denoising using 1D dictionary learning 
% Prepared By Yangkang Chen and Hang Wang
% Dec, 09, 2018
%
% Key Reference
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51-KS61.
% 
% NOTE: that this is an erratic noise example
% The algorithm is described in 
% Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iteratively robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51-KS61.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2020, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 222, 1846?1863. 
% Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2020, Erratic-noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, in press
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

n1=256;
r=zeros(256,1);
r(100)=1;
r(180)=0.7;
% r(240)=0.3;
w=yc_ricker(20,0.004,0.1);
d=conv(r,w,'same');
d=yc_phase_correction(d,90);
figure;plot(d);

d0=d;
randn('state',2015161718);
d=d0+0.1*randn(size(d0));

t=(1:256)*2;
 figure('units','normalized','Position',[0.2 0.4 0.4, 0.4]);
subplot(2,1,1);plot(t,d0,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
subplot(2,1,2);plot(t,d,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);

randn('state',200000124);
err_n=100*randn(40,1);
d_err=d;d_err(121:160)=d_err(115:154)+err_n;

figure('units','normalized','Position',[0.2 0.4 0.4, 0.4]);
subplot(2,1,1);plot(t,d0,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
subplot(2,1,2);plot(t,d_err,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);

% tmp=mfcyk(d_err,20,1,1);
% tmp=mfcyk(tmp,20,1,1);
tmp=yc_mfs(d_err,20,1,0,3);
figure;plot(tmp);

tmp=err_n;
randn('state',232043124);
tmp(find(err_n>0.2 | err_n<-0.2 ))=0.1*randn(size(find(err_n>0.2| err_n<-0.2)));
d_pre=d;d_pre(121:160)=tmp;
figure;plot(d_pre);
dd=d_pre;
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

niter=10;
m1=0;
m2=0;
for it=1:niter

n1=d_err-m1;
n2=d_err-m2;

% randn('state',232043124);
n1(find(n1>1.00 | n1<-1.00 ))=0.1*randn(size(find(n1>1.00| n1<-1.00)));
% randn('state',2320431245);
n2(find(n2>1.00 | n2<-1.00 ))=0.1*randn(size(find(n2>1.00| n2<-1.00)));

p1=m1+0.5*n1;
p2=m2+0.5*n2;

dd1=p1;
dd2=p2;
%% decompose the image into patches:
X=yc_patch(dd1,1,16,1,8,1);
XX=yc_patch(dd2,1,16,1,8,1);

%% for ppt
% figure;plot(70:150,d(70:150),'r','linewidth',4);
% Dt=X(:,3:14);Dt(:,2:2:end)=zeros(size(Dt(:,2:2:end)));
% figure;yc_wigbh(Dt);
% figure('units','normalized','Position',[0.2 0.4 0.3, 0.8],'color','w');wigbh_tmp(X(:,9:14)+0.25*randn(16,6),1:6,1:16,1.5);
% 
% ylim([0,7]);

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




% figure('units','normalized');
% imagesc(G);colormap(jet);caxis([-1,2.0]);colorbar;
% ylabel('Atom NO','Fontsize',16);
% xlabel('Patch NO','Fontsize',16);
% title('Coefficients Matrix','Fontsize',16);
% set(gca,'Linewidth',1.5,'Fontsize',16);
% print(gcf,'-depsc','-r300','syn_G.eps');

%% K-SVD

params.data = XX;
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

m1=d1;
m2=d2;

fprintf("iteration = %d/%d,snr(m1)=%g, snr(m2)=%g, \n",it,niter,yc_snr(d0,d1),yc_snr(d0,d2));

end



%%plot

dt=2;t=[0:n1-1]*dt;
% yc_snr(d0,d2)   %7.48

figure;%imagesc([d,d1,d-d1]);
subplot(3,1,1);plot(d);ylim([-1,1]);
subplot(3,1,2);plot(d2);ylim([-1,1]);
subplot(3,1,3);plot(d-d2);ylim([-1,1]);


%  figure('units','normalized','Position',[0.2 0.4 0.4, 0.4]);
% subplot(2,1,1);plot(t,d0,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
% ylabel('Amplitude','Fontsize',16);
% set(gca,'Linewidth',1.5,'Fontsize',16);
% title('(a)','Fontsize',16);
% 
% subplot(2,1,2);plot(t,d,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
% ylabel('Amplitude','Fontsize',16);
% xlabel('Time (ms)','Fontsize',16);
% set(gca,'Linewidth',1.5,'Fontsize',16);
% title('(b)','Fontsize',16);
% print(gcf,'-depsc','-r300','syn_d.eps');

%BP
d3=yc_bp(d,0.002,0,0.5,100,105);

%EMD
d4=yc_micro_emd(d);  

yc_snr(d0,d2)    %6.6491
yc_snr(d0,d1)   %10.75
yc_snr(d0,d3)   %7.5490
yc_snr(d0,d4)   %4.07
yc_snr(d0,d)   %3.25

figure('units','normalized','Position',[0.2 0.4 0.4, 0.6]);
subplot(4,1,1);plot(t,d0,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(4,1,2);plot(t,d_err,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);

subplot(4,1,3);plot(t,d2,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);

subplot(4,1,4);plot(t,d1,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(d)','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
print(gcf,'-depsc','-r300','syn_dn_err.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.6]);
subplot(4,1,1);plot(t,d2,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(a)','Fontsize',16);

subplot(4,1,2);plot(t,d1,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(b)','Fontsize',16);

subplot(4,1,3);plot(t,d_err-d2,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(c)','Fontsize',16);

subplot(4,1,4);plot(t,d_err-d1,'b','linewidth',2);ylim([-1,1]);xlim([0,t(end)]);
ylabel('Amplitude','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
title('(d)','Fontsize',16);
xlabel('Time (ms)','Fontsize',16);
print(gcf,'-depsc','-r300','syn_dn_err_n.eps');


