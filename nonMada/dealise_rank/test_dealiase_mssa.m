%% This is a demo for de‐aliased and de‐noised Cadzow filtering for seismic data reconstruction
% Prepared by Weilin Huang and Yangkang Chen
% June, 2021
% 
% Reference
% Huang, W., D. Feng, and Y. Chen, 2020, De‐aliased and de‐noise Cadzow
% filtering for seismic data reconstruction, Geophysical Prospecting, 68,
% 443-571.

clear;close all;clc;

%% create data
wt=300;
wx=80;
p=2.2;
a1=zeros(wt,wx);
[n,m]=size(a1);
a2=a1;
a3=a1;
a4=a1;
a5=a1;
a6=a1;
a7=a1;
a8=a1;
a9=a1;
a10=a1;
a11=a1;
a12=a1;
a13=a1;
a14=a1;
k=0;
a=0.1;
b=1;
for t=-0.054:0.001:0.056
    k=k+1;
    b1(k)=(1-2*(pi*35*t).^2).*exp(-(pi*35*t).^2);
    b2(k)=(1-2*(pi*20*t).^2).*exp(-(pi*20*t).^2);
    b3(k)=(1-2*(pi*60*t).^2).*exp(-(pi*60*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
randn('state',1);
c=0.2*round(randn(1,200));
c=0.2*ones(1,200);
c([1:10:199])=1;
c([2:20:200])=-1;
%bnew=conv(c,b1);
bnew=sin(0:0.3:20*pi);
kk=size(bnew,2);
for j=1:20
for i=1:m
%     t1(i)=round(a*i^2+b);
%     t2(i)=round(0.5*i^2+100);
%     t3(i)=round(0.02*i^2+200);
%     t4(i)=round(5*sin(i/3)+350);
%     a1(t1(i):t1(i)+k-1,i)=b1;
%     a2(20+i:20+k-1+i,i)=b2;
%     a3(40+2*i:40+k-1+2*i,i)=b3;
%     a4(40+2*i:40+k-1+2*i,i)=b4;
%     a5(300-i:300+k-1-i,i)=b2;
%     a6(470-3*i:470+k-1-3*i,i)=b3;
%     a7(t2(i):t2(i)+k-1,i)=b3;
%     a8(t3(i):t3(i)+k-1,i)=b4;
%     a9(30:30+k-1,i)=b4;
%  a14(t4(i):t4(i)+k-1,i)=b2;

%   t1(i)=round(0.005*(i-1)^2+110);
%   t11(i)=1;
%   t2(i)=round(0.009*(i-1)^2+50);
%   t3(i)=round(0.003*(i-1)^2+170);
  t4(i)=round(-1*i+185-1*j);
  t5(i)=round(2*i+15);
  t6(i)=round(50+j*3);
%   a1(t1(i):t1(i)+k-1,i)=b1;
%   % a1(t11(i):t11(i)+k-1,i)=b1;
%   a2(t2(i):t2(i)+k-1,i)=b1; 
%   a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b2;
  a5(t5(i):t5(i)+k-1,i)=-b2;
  a6(t6(i):t6(i)+k-1,i)=b1;
end

%% add band-limitted noise and decimate data
randn('state',201315+j);
%shot=0*randn(wt,wx)+a2(1:wt,:)+a1(1:wt,:)+a3(1:wt,:)+a4(1:wt,:);%+a5(1:500,:)+a6(1:500,:)+a7(1:500,:)+a8(1:500,:)+a9(1:500,:)+a10(1:500,:)+a11(1:500,:)+a12(1:500,:)+a13(1:500,:)+a14(1:500,:);%+randn(100,20,10);
noise=1*randn(wt,wx);
% f_n=fft(noise);
% f_n(1:3,:)=0;
% f_n(15:100,:)=0;
% noise=real(ifft(f_n));
[o] =  yc_bp(noise,0.001,2,10,90,100);
shotn(:,:,j)=0.6*o+a5(1:wt,:)+a4(1:wt,:)+a6(1:wt,:);
shotc(:,:,j)=a5(1:wt,:)+a4(1:wt,:)+a6(1:wt,:);
end
shotn2=shotn;
jj=1;
for j=2:2:20
    shotn(:,:,j)=0;
    ii=1;
for i=2:2:wx
    shotn(:,i,j-1)=0;
    shotn_in(:,ii,jj)=shotn(:,i-1,j-1);
    ii=ii+1;
end
jj=jj+1;
end

% load('clean_MAP')
% hwl_3D_image(shotn_in,2,'Clean data','x','y','Time (ms)',MAP_clean_model,d,u)
% hwl_3D_image(shotn,2,'Clean data','x','y','Time (ms)',MAP_clean_model,d,u)
% hwl_3D_image(shotc,2,'Clean data','x','y','Time (ms)',MAP_clean_model,d,u)

% load('clean_MAP')
% hwl_3D_image(shotc,2,'Clean data','x','y','Time (ms)',MAP_clean_model,d,u)
% set(gcf,'Position',[100 100 880 850]);
% v3d_slice(4/5,4/5,4/5)
% hwl_3D_image(shotn,2,'Clean data','x','y','Time (ms)',MAP_clean_model,d,u)
% set(gcf,'Position',[100 100 880 850]);
% v3d_slice(4/5,4/5,4/5)
% y=m;
% z=20;
% x=wt;
% rand('state',199079);
% yy=randperm(y*z-1);
% % for i=1:round(y*z*0.3)
% %     cy=rem(yy(i),y);
% %     cz=(yy(i)-cy)/y;
% %     shot(:,cy+1,cz+1)=0;
% % end
% shot_r=reshape(shot,[x,y*z]);
% for i=1:2:y*z
%     shot_r(:,i:i)=zeros(x,1);
% end
% shot=reshape(shot_r(:,1:y*z),[x,y,z]);

%% Doing the reconstruction
d0=shotn;
d1=hwl_double_MSSA_noise(d0,1,50,0.004,3,2,10,2);
dn=shotn2;
figure;imagesc(d1(:,:,1));
figure;imagesc(d0(:,:,1));

d2=hwl_double_MSSA_noise(dn,1,50,0.004,3,2,10,2);%takes about 20 minutes
figure;imagesc(dn(:,:,1));
figure;imagesc(d2(:,:,1));

