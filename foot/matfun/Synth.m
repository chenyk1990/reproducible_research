function Synth(clean,nft,noisy,d_sgrdl,d_mf,d_dl,d_drr)
% Author      : Wei Chen and Yangkang Chen
%               Yangtze University and Zhejiang University
%         
% Date        : Feb, 2021
%
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%
%  Copyright (C) 2021 Yangtze University and Zhejiang University
%  Copyright (C) 2021 Wei Chen and Yangkang Chen
%  
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
%  
%  Reference:
%  Chen et al., 2021, Statistics-guided residual dictionary learning for footprint noise removal, IEEE TGRS, doi: 10.1109/TGRS.2021.3070903.
%  Zhou et al., 2021, Statistics-guided dictionary learning for automatic coherent noise suppression, IEEE TGRS, doi: 10.1109/TGRS.2020.3039738.
%  Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE TGRS, doi: 10.1109/TGRS.2020.3030740.
%  Chen, 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, GJI, 222, 1717-1727.
% 
% Synthetic data for footprint noise removal
% clc;clear;close all;

%% generate synthetic data
%plane1
n1=300;n2=40;n3=40;
t1=150;p11=0;p12=0;
plane1=zeros(n1,n2,n3);
for i2=1:n2
    for i3=1:n3
        
        t=floor(t1+p11*i2+p12*i3);
        plane1(t,i2,i3)=1;
    end
end
plane1=ricker_mada(plane1,0.002,10);

%plane2
n1=300;n2=40;n3=40;
t2=50;p21=2;p22=2;
plane2=zeros(n1,n2,n3);
for i2=1:n2
    for i3=1:n3
        
        t=floor(t2+p21*i2+p22*i3);
        plane2(t,i2,i3)=1;
    end
end
plane2=ricker_mada(plane2,0.002,10);

%%plane3
n1=300;n2=40;n3=40;
t3=230;p31=-2;p32=-2;
plane3=zeros(n1,n2,n3);
for i2=1:n2
    for i3=1:n3
        
        t=floor(t3+p31*i2+p32*i3);
        plane3(t,i2,i3)=1;
    end
end
plane3=ricker_mada(plane3,0.002,10);

dc=plane1+plane2+plane3;
dc=yc_scale(dc,3);
%figure;imagesc([dc(:,:)]);colormap(seis);
%figure;imagesc([dc(:,:,5)]);colormap(seis);

%% adding random noise
randn('state',201314);
var=0.2;
% dn=d+var*randn(size(d));
n=0.02*randn(size(dc));

%% adding footprint noise
[nt,nx,ny]=size(dc);
x=[0:nx-1];
y=sin(10*x);
%figure;plot(y);
noise2d=y(:)*ones(1,ny);
%
%figure;imagesc(noise2d);

noise3d=zeros(nt,nx,ny);
for it=0:nt-1
    
    noise3d(it+1,:,:) = noise2d*0.1^(it/(nt-1));
end

dn=dc+noise3d*0.1+n;

%% see the data
%figure;imagesc([squeeze(dc(100,:,:)),squeeze(dn(100,:,:))]);colormap(seis);
%figure;imagesc([dc(:,:,9),dn(:,:,9)]);colormap(seis);

%% denoise by DRR
flow=0;fhigh=125;dt=0.004;N=3;verb=1;
% d1=fxydmssa(dn(:,:,:),flow,fhigh,dt,N,6,verb);
% figure;imagesc([dc(:,:,9),dn(:,:,9),d1(:,:,9),dn(:,:,9)-d1(:,:,9)]);

%% denoise by KSVD
d=squeeze(dn(100,:,:));

% d=yc_transp(dn,13);

l1=8;l2=8;s1=4;s2=4;
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

%% Denoising by KSVD
%param naming following Chen, 2017, GJI; Zhou et al., 2020
K=64;
param.T=3;      %sparsity level
param.D=DCT;    %initial D
param.niter=10; %number of K-SVD iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
%param.exact:   Exact K-SVD update or approximate
param.K=64;     %number of atoms, dictionary size
%for X=DG
%size of X: MxN
%size of D: MxK
%size of G: KxN

%% Option 1: denoise only using the integrated function
param=struct('T',3,'niter',10,'mode',1,'K',64,'D',DCT);
mode=1;l1=8;l2=8;s1=4;s2=4;perc=7;
% d1=yc_ksvd_denoise(d,mode,[l1,l2,1],[s1,s2,1],perc,param);
% figure;imagesc([dc,d,d1,d-d1]);colormap(seis);
% yc_snr(dc,d1)
% 
% %% Option 2: without initialized param.D
% param=rmfield(param,'D');%param=struct('T',3,'niter',10,'mode',1,'K',64);
% mode=1;l1=8;l2=8;s1=4;s2=4;perc=7;
% d1=yc_ksvd_denoise(d,mode,[l1,l2,1],[s1,s2,1],perc,param);
% figure;imagesc([dc,d,d1,d-d1]);colormap(seis);
% yc_snr(dc,d1)
% 
%% compare performance of two dictionaries
T=3;perc=100;
[n1,n2]=size(d);
%% residual learning
X=yc_patch(d,mode,l1,l2,s1,s2);
[D,G]=yc_ksvd(X,param);
G2=yc_ompN(D,X,T);
% G2=yc_pthresh(G2,'ph',perc);
X2=D*G2;
d3=yc_patch_inv(X2,mode,n1,n2,l1,l2,s1,s2);
% d3=yc_transp(d3,13);
% figure;imagesc([squeeze(dn(100,:,:)),squeeze(d3(100,:,:))]);colormap(seis);
%figure;imagesc([d,d3,d-d3]);colormap(seis);

%% filtered dictionary
%stronger
DD=zeros(size(D));
for ia=1:64
    DD(:,ia)=reshape(yc_mfs(reshape(D(:,ia),l1,l2),2,1,1,4),l1*l2,1);
    DD(:,ia)=reshape(yc_mfs(reshape(DD(:,ia),l1,l2),2,1,2,4),l1*l2,1);
end
%weaker
DDD=zeros(size(D));
for ia=1:64
    DDD(:,ia)=reshape(yc_mfs(reshape(D(:,ia),l1,l2),2,1,1,2),l1*l2,1);
    DDD(:,ia)=reshape(yc_mfs(reshape(DDD(:,ia),l1,l2),2,1,2,2),l1*l2,1);
end

%% residual DL
G22=yc_ompN([DD,D-DD],X,3);
X22=DD*G22(1:K,:);
d33=yc_patch_inv(X22,mode,n1,n2,l1,l2,s1,s2);
%figure;imagesc([squeeze(dc(100,:,:)),d,d33,d-d33]);colormap(seis);
yc_snr(squeeze(dc(100,:,:)),squeeze(d33(:,:)))

%% residual DL + statistics guided
natom=K;
% k=kurtosis(D);%16.6513 N=21; 16.9889 N=22; 17.1545 N=23; 17.3897 N=24; 17.6918 N=25; 17.1196 N=26
% k=yc_kurtosis2(D,l1,l2);
k=yc_var2(D,l1,l2);
[ks,ii]=sort(k,'descend');
% figure;stem(ks);
% ks_ratio=[ks(1:end-1)./ks(2:end)];% ks_dif=diff(ks);figure;stem(ks_dif);
% figure;stem(ks_ratio);
perc=15;%works fine
% perc=90;%works perfectly? I think so.
tt=round(natom*(100-perc)/100);
inds=ii(1:tt);
% 
% D_o1=D;
% D_o1(:,inds)=0;
% figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
% for ia=1:64
%     subplot(8,8,ia);imagesc(reshape(D(:,ia),l1,l2));colormap(jet);
%     set(gca,'Linewidth',1.5,'Fontsize',16);
%     set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
% end
% 
% figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
% for ia=1:64
%     subplot(8,8,ia);imagesc(reshape(D_o1(:,ia),l1,l2));colormap(jet);
%     set(gca,'Linewidth',1.5,'Fontsize',16);
%     set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
% end

% DD_o=D;
DD_o=DDD;
DD_o(:,inds)=DD(:,inds);
G22_o=yc_ompN([DD_o,D-DD_o],X,3);
X22_o=DD_o*G22_o(1:K,:);
d33_o=yc_patch_inv(X22_o,mode,n1,n2,l1,l2,s1,s2);
%figure;imagesc([squeeze(dc(100,:,:)),d,d33_o,d-d33_o]);colormap(seis);
yc_snr(squeeze(dc(100,:,:)),squeeze(d33_o(:,:)))


%% without residual DL
G222=yc_ompN([DD],X,12);
X222=D*G222(1:K,:);
d333=yc_patch_inv(X222,mode,n1,n2,l1,l2,s1,s2);
%figure;imagesc([d,d333,d-d333]);colormap(seis);

%with mf
d44=yc_mfs(d,4,1,1,1);
%figure;imagesc([[d,d33,d-d33];[d,d44,d-d44]]);colormap(gray);colormap(seis);
yc_snr(squeeze(dc(100,:,:)),squeeze(d44(:,:)))

%% plot the first 64 atoms
%figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
%for ia=1:64
%    subplot(8,8,ia);imagesc(reshape(DCT(:,ia),l1,l2));colormap(jet);
%    set(gca,'Linewidth',1.5,'Fontsize',16);
%    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
%end
% print(gcf,'-depsc','-r300','syn_atoms.eps');

%figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
%for ia=1:64
%    subplot(8,8,ia);imagesc(reshape(D(:,ia),l1,l2));colormap(jet);
%    set(gca,'Linewidth',1.5,'Fontsize',16);
%    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
%end

%figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
%for ia=1:64
%    subplot(8,8,ia);imagesc(reshape(DD(:,ia),l1,l2));colormap(jet);
%    set(gca,'Linewidth',1.5,'Fontsize',16);
%    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
%end

%figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
%for ia=1:64
%    subplot(8,8,ia);imagesc(reshape(D(:,ia)-DD(:,ia),l1,l2));colormap(jet);
%    set(gca,'Linewidth',1.5,'Fontsize',16);
%    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
%end

%% with statistics
%figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
%for ia=1:64
%    subplot(8,8,ia);imagesc(reshape(DD_o(:,ia),l1,l2));colormap(jet);
%    set(gca,'Linewidth',1.5,'Fontsize',16);
%    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
%end

%figure('units','normalized','Position',[0.2 0.4 0.6, 0.8]);
%for ia=1:64
%    subplot(8,8,ia);imagesc(reshape(D(:,ia)-DD_o(:,ia),l1,l2));colormap(jet);
%    set(gca,'Linewidth',1.5,'Fontsize',16);
%    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
%end


%% for all time slices
d2=zeros(size(dc)); 	%sgrdl
d22=zeros(size(dc));	%mf
d222=zeros(size(dc));	%dl
d2222=zeros(size(dc));	%mssa

for i1=1:nt
diter=squeeze(dn(i1,:,:));
Xiter=yc_patch(diter,mode,l1,l2,s1,s2);
if i1<=30
G22iter=yc_ompN([DD_o,D-DD_o],Xiter,2); 
else
G22iter=yc_ompN([DD_o,D-DD_o],Xiter,6);
end
X22iter=DD_o*G22iter(1:K,:);
d33iter=yc_patch_inv(X22iter,mode,n1,n2,l1,l2,s1,s2);
% figure;imagesc([diter,d33iter,diter-d33iter]);colormap(seis);
d2(i1,:,:)=d33iter;
d22(i1,:,:)=yc_mfs(diter,4,1,1,1);

G222iter=yc_ompN(D,Xiter,6); 
X222iter=D*G222iter;
d222(i1,:,:)=yc_patch_inv(X222iter,mode,n1,n2,l1,l2,s1,s2);
fprintf('i1=%d/%d is done\n',i1,nt);
end

%% denoise by DRR
flow=0;fhigh=125;dt=0.004;N=3;verb=1;
d2222=fxydmssa(dn,flow,fhigh,dt,N,6,verb);

% figure;imagesc([dc(:,:,9),dn(:,:,9),d1(:,:,9),dn(:,:,9)-d1(:,:,9)]);colormap(seis);yc_snr(dc,d1,3)
%figure;imagesc([dc(:,:,9),dn(:,:,9),d2(:,:,9),dn(:,:,9)-d2(:,:,9)]);colormap(seis);yc_snr(dc,d2,3)
%21.33 (without statistics
%figure;imagesc([dc(:,:);dn(:,:);d2(:,:)]);colormap(seis);

% for the first few
% diter=squeeze(dn(10,:,:));
% Xiter=yc_patch(diter,mode,l1,l2,s1,s2);
% G22iter=yc_ompN([DD,D-DD],Xiter,1);
% X22iter=DD*G22iter(1:K,:);
% d33iter=yc_patch_inv(X22iter,mode,n1,n2,l1,l2,s1,s2);
% figure;imagesc([diter,d33iter,diter-d33iter]);colormap(seis);



rsf_create(clean,size(dc)');
rsf_write(dc,clean);

rsf_create(nft,size(noise3d*0.1)');
rsf_write(noise3d*0.1,nft);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(d_sgrdl,size(d2)');
rsf_write(d2,d_sgrdl);

rsf_create(d_mf,size(d22)');
rsf_write(d22,d_mf);

rsf_create(d_dl,size(d222)');
rsf_write(d222,d_dl);

rsf_create(d_drr,size(d2222)');
rsf_write(d2222,d_drr);

























