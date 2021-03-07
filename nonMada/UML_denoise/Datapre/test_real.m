% clear all
% clc
% Step1: Pre-processing
% Step2: Tensorflow
% Step3: Post-processing

%% Pre-processing
% Dreal=readdat('real3d-10.dat',512,128); % in the paper
Dreal=readdat('real3d-32.dat',512,128); % a slightly different data

% Normalization
d1=Dreal;
min1=min(min(min(d1)));
max1=max(max(max(d1)));
d2=(d1-min1)./(max1-min1);

% Randomly sampled patch
AA=patchmi(d2,10000,40,40);
AAA=patchmi(d2,10000,40,40);

% Regularly sampled patch
A0=patch_cyk(d2,1,40,40,20,20)';


d2_tmp=zeros(520,140);
d2_tmp(1:512,1:128)=Dreal;
figure;imagesc(d2_tmp);

figure;imagesc(A0');
AA(1:150,:)=A0;

% Save data for tensorflow
% Train dataset
writdat1('Real2-10000-1600-n-edge.dat',AA,10000,1600)
writdat1('Real2-10000-1600-n-edge2.dat',AAA,10000,1600)

% Test dataset
writdat1('Real2-150-1600-edge.dat',A0,150,1600)

% %% Post-processing
dresaultreal=readdat1('denoised-real2-150-1600-edge.dat',150,1600);
DresultReal=yc_patch_inv(dresaultreal',1,512,128,40,40,20,20);

dresaultreal2=readdat1('denoised-real2-150-1600-edge2.dat',150,1600);
DresultReal2=yc_patch_inv(dresaultreal2',1,512,128,40,40,20,20);


figure;imagesc(A0');colormap(seis);caxis([0,1]);
% figure('units','normalized','Position',[0.2 0.4 0.55, 0.35],'color','w');
ylabel('Pixel NO','Fontsize',20);
xlabel('Sample NO','Fontsize',20);
set(gca,'Linewidth',1.5,'Fontsize',20);
print(gcf,'-depsc','-r400','real_A0.eps');

figure;imagesc(dresaultreal');colormap(seis);caxis([0,1]);
ylabel('Pixel NO','Fontsize',20);
xlabel('Sample NO','Fontsize',20);
set(gca,'Linewidth',1.5,'Fontsize',20);
print(gcf,'-depsc','-r400','real_A1.eps');

figure;imagesc(dresaultreal2');colormap(seis);caxis([0,1]);
ylabel('Pixel NO','Fontsize',20);
xlabel('Sample NO','Fontsize',20);
set(gca,'Linewidth',1.5,'Fontsize',20);
print(gcf,'-depsc','-r400','real_A2.eps');


% Inverse Normalization
d3=(max1-min1)*DresultReal+min1;
d33=(max1-min1)*DresultReal2+min1;

% Save final result
writdat('Real2-512-128-175-1600-denoised-SAE-edge.dat',d3,512,128)

figure;imagesc([d1,d3,d1-d3]);colormap(seis);
figure;imagesc([d1,d33,d1-d33]);colormap(seis);


x=[1:128];
t=[0:511]*0.004;

x0=[1:140];
t0=[0:520]*0.004;

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8]);
imagesc(x0,t0,d2_tmp);colormap(seis);caxis([-3,3]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
%title('Noise','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
print(gcf,'-depsc','-r200','real_edge0.eps');

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8]);
imagesc(x,t,d3);colormap(seis);caxis([-3,3]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
%title('Noise','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
print(gcf,'-depsc','-r200','real_edge1.eps');

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8]);
imagesc(x,t,d33);colormap(seis);caxis([-3,3]);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Trace','Fontsize',20,'fontweight','bold');
%title('Noise','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
print(gcf,'-depsc','-r200','real_edge2.eps');

% %
% %
% % figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
% % imagesc(d1);colormap(seis);caxis([-3,3]);
% % ylabel('Time (s)','Fontsize',20);
% % xlabel('Trace','Fontsize',20);
% % title('Real data','Fontsize',20);
% % set(gca,'Linewidth',2,'Fontsize',20);
% % print(gcf,'-depsc','-r200','real_dn.eps');
% %
% % figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
% % imagesc(d3);colormap(seis);caxis([-3,3]);
% % ylabel('Time (s)','Fontsize',20);
% % xlabel('Trace','Fontsize',20);
% % title('Denoised data','Fontsize',20);
% % set(gca,'Linewidth',2,'Fontsize',20);
% % print(gcf,'-depsc','-r200','real_dsae.eps');
% %
% % figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
% % imagesc(d1-d3);colormap(seis);caxis([-3,3]);
% % ylabel('Time (s)','Fontsize',20);
% % xlabel('Trace','Fontsize',20);
% % title('Noise','Fontsize',20);
% % set(gca,'Linewidth',2,'Fontsize',20);
% % print(gcf,'-depsc','-r200','real_dsae_n.eps');
% %
% %
% % %% compared with FX SSA
% % d4=fxymssa(d1,0,100,0.004,5,1);
% % figure;imagesc([d1,d4,d1-d4]);colormap(seis);
% %
% % figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
% % imagesc(d4);colormap(seis);caxis([-3,3]);
% % ylabel('Time (s)','Fontsize',20);
% % xlabel('Trace','Fontsize',20);
% % title('Denoised data','Fontsize',20);
% % set(gca,'Linewidth',2,'Fontsize',20);
% % print(gcf,'-depsc','-r200','real_dssa.eps');
% %
% % figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
% % imagesc(d1-d4);colormap(seis);caxis([-3,3]);
% % ylabel('Time (s)','Fontsize',20);
% % xlabel('Trace','Fontsize',20);
% % title('Noise','Fontsize',20);
% % set(gca,'Linewidth',2,'Fontsize',20);
% % print(gcf,'-depsc','-r200','real_dssa_n.eps');
%
%
%
% %% MSSA DRR
% % N3=6;
% % N4=7;
% % K2=5;
% % flow=0;fhigh=250;dt=0.001;verb=1;
% %
% % Dmssa2=fxymssa(d2,flow,fhigh,dt,N3,verb);
% % Ddmssa2=fxydmssa(d2(:,:,:),flow,fhigh,dt,N4,K2,verb);
% %
% % Ddmssa2=(max1-min1)*Ddmssa+min1;
% % Dmssa2=(max1-min1)*Dmssa+min1;
% %
% % writdat('Real-512-128-dmssa-R-edge.dat',Ddmssa2,512,128)
% % writdat('Real-512-128-mssa-R-edge.dat',Dmssa2,512,128)
