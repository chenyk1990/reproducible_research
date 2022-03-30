clc;clear;close all;


addpath(genpath('~/chenyk.data/pywinml/matfun'));


%%
%% create data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
    t1(i)=round(140);
    t3(i)=round(-6*i+180);
    t4(i)=round(6*i+10);
    a1(t1(i):t1(i)+k-1,i)=b1;
    a3(t3(i):t3(i)+k-1,i)=b1;
    a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
        t4(i)=round(6*i+10+3*j);
        a4(t4(i):t4(i)+k-1,i)=b1;
        
        t1(i)=round(140-2*j);
        a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
dc=yc_scale(plane3d,3);
dc=dc(51:225,:,:);
figure;imagesc(reshape(dc,175,20*20));colormap(seis);

%% adding noise
randn('state',201617);
dn=dc+0.02*randn(size(dc));
[n1,n2,n3]=size(dn);
figure;imagesc([dc(:,:,10),dn(:,:,10)]);colormap(seis);

%% SGK
l1=16;l2=4;l3=4;c1=16;c2=4;c3=4;
s1=1;s2=1;s3=1;

% l1=8;l2=8;l3=8;c1=l1;c2=l2;c3=l3;
% s1=l1/2;s2=l2/2;s3=l3/2;
perc=1;
[DCT]=yc_initD([l1,l2,l3],[c1,c2,c3]);
param.mode=1;   %1: sparsity; 0: error
param.niter=10; %number of K-SVD iterations to perform; default: 10
param.D=DCT;    %initial D
param.T=3;      %sparsity level

X=yc_patch3d(dn,1,l1,l2,l3,s1,s2,s3);
size(X)
% [Dsgk,Gsgk]=yc_sgk(X,param);
% Gsgkc=Gsgk;
% Gsgk=yc_pthresh(Gsgkc,'ph',perc);
% X1=Dsgk*Gsgk;
% d1=yc_patch3d_inv(X1,1,n1,n2,n3,l1,l2,l3,s1,s2,s3);
% figure;imagesc([dc(:,:,10),dn(:,:,10),d1(:,:,10)]);colormap(seis);
% yc_snr(dc,dn,2)
% yc_snr(dc,d1,2)
%% MLP
%size of X: Patch_size x N_patches
%size of x,y,y1: N_patches x Patch_size
x=X';
y=X';

% figure(1);
% scatter(X(:,1),X(:,2),20,y);
% hold on;

% A NETWORK WITH A HIDDEN LAYER OF SIZE 3
nn_hdim=128;niter=2000;verb=1;lrate=0.01;nepoch=200;
model = nns_build_model(x,y,nn_hdim,niter,lrate,nepoch,verb);
y1=nns_predict(model,x); %size of x,y1: N_patches x Patch_size
d2=yc_patch3d_inv(y1',1,n1,n2,n3,l1,l2,l3,s1,s2,s3);
figure;imagesc([dc(:,:,10),dn(:,:,10),d2(:,:,10)]);colormap(seis);
yc_snr(dc,d2,2)

