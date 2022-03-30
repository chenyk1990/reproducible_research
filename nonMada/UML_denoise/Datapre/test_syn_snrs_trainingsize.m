clc;clear;close all;

d=readdat('SaltShot-5m-1ms-1500-180.dat',1500,180);

p1=round(size(d,1)*0.37)-3;
data1=d(1600-p1+1:1500,61:180); % Single event example
dc=d(101:p1,61:180); % Complex example

% d1=data1;
% min1=min(min(min(d1)));
% max1=max(max(max(d1)));
% d2=(d1-min1)./(max1-min1);

randn('state',201809);
n=randn(size(dc));
dn=dc+0.1*n;

min1=min(min(min(dn)));
max1=max(max(max(dn)));
dn_tmp=(dn-min1)./(max1-min1);

Ntrain=1000;
Ntrain=6000;
s1=40;
% s2=round(s1/2);
s2=30;

A=patchmi(dn_tmp,Ntrain,s1,s1);
A0=yc_patch(dn_tmp,1,s1,s1,s2,s2)';
nA0=size(A0,1)
A(1:nA0,:)=A0;

writdat1('syn_A.dat',A,Ntrain,s1*s1)
writdat1('syn_A0.dat',A0,nA0,s1*s1)


%% the following is run after training and denoising
tmp=readdat1('syn_A1.dat',nA0,s1*s1);
dtmp=yc_patch_inv( tmp',1,452,120,s1,s1,s2,s2);

d1=(max1-min1)*dtmp+min1;

figure;imagesc([dn,d1,dn-d1]);
figure;yc_imagesc([dn,d1,dn-d1]);

yc_snr(dc,dn)
yc_snr(dc,d1)


%%
N=[1000:1000:6000];
snrs=[7.56,10.46,10.91,11.98,12.54,12.55];

figure;
plot(N,snrs,'or-','linewidth',2);
ylabel('SNR (dB)','Fontsize',16,'fontweight','bold');
xlabel('Training Data Size','Fontsize',16,'fontweight','bold');
title('SNR V.S. Training Data Size','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
print(gcf,'-depsc','-r400','syn_size.eps');

%%
S=[20,30,40,50,60];
snrs2=[0.43,4.16,10.38,6.23,3.72];

figure;
plot(S,snrs2,'or-','linewidth',2);
ylabel('SNR (dB)','Fontsize',16,'fontweight','bold');
xlabel('Patch Size','Fontsize',16,'fontweight','bold');
title('SNR V.S. Patch Size','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
print(gcf,'-depsc','-r400','syn_patch_size.eps');

%%
S2=[2,4,10,20,30];
snrs3=[13.95,13.63,11.03,10.38,4.16];

figure;
plot(S2,snrs3,'or-','linewidth',2);ylim([0,15]);xlim([1,30]);
ylabel('SNR (dB)','Fontsize',16,'fontweight','bold');
xlabel('Shift Size','Fontsize',16,'fontweight','bold');
title('SNR V.S. Shift Size','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
print(gcf,'-depsc','-r400','syn_shift_size.eps');

%N_train
%7.56 dB, 1000
%10.46 dB, 2000
%10.91 dB, 3000
%11.98 dB, 4000
%12.54 dB, 5000
%12.55 dB, 6000
%12.47 dB, 7000
%12.30 dB, 8000
% dn: -5.05 dB

%Size Patch
%0.43 dB, 20, 6000
%4.16 dB, 30, 6000
%10.38 dB, 40, 6000
%6.23 dB, 50, 6000
%3.72 dB, 60, 6000

%Overlapping size
%4.16 dB, 40, 6000, 30(shift)
%10.38 dB, 40, 6000, 20(shift)
%11.03, 40, 6000, 10(shift)
%13.63 , 40, 6000, 4(shift)
%13.95 dB, 40, 6000, 2(shift)




