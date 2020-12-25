function linear(clean,noisy,decimated,lmafit,svd,lanb,rsvd)
%% data preparition and parameterization 
addpath(genpath('../matfun'));
load ../../drr5d/matfun/yc_synth5d.mat
[nt,nhx,nhy,nx,ny]=size(data5d);
d=data5d;

dt=0.004;
nf=2^nextpow2(nt);
flow=0;
fhigh=125;
N1=3; % rank for svd and lmafit
N2=5; % rank for lanb
N3=3; % rank for rsvd
P=6; % trade-off for rsvd
K=3; % damped factor for dmssa
eps=0.00001;
Niter=10;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing

% figure;
% imagesc(squeeze(d(:,5,5,:,6)));

% d=permute(d,[1 4 5 2 3]); % switch hx,hy with x,y

%% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));

%% decimating data
ratio=0.3;
mask=yc_genmask(reshape(dn,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
dn0=dn.*mask;

%% data processing - LR5D
tic
[d_recon_lmafit] = lr5d_lmafit(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N1,K,eps,Niter,a,mask);
toc
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_lmafit,nt*nhx*nhy,nx,ny),2)

tic
[d_recon_svd] = lr5d_svd(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N1,K,eps,Niter,a,mask);
toc
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_svd,nt*nhx*nhy,nx,ny),2)

tic
[d_recon_lanb] = lr5d_lanb(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N2,K,eps,Niter,a,mask);
toc
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_lanb,nt*nhx*nhy,nx,ny),2)

tic
[d_recon_rsvd] = lr5d_rsvd(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N3,K,P,eps,Niter,a,mask);
toc
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
yc_snr(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_rsvd,nt*nhx*nhy,nx,ny),2)

% figure;
% imagesc([squeeze(d(:,5,5,:,6)),squeeze(dn(:,5,5,:,6)),squeeze(d_recon_lmafit(:,5,5,:,6)),...
%     squeeze(d_recon_lanb(:,5,5,:,6)),squeeze(d_recon_rsvd(:,5,5,:,6)),...
%     squeeze(d_recon_lmafit(:,5,5,:,6))-squeeze(d(:,5,5,:,6)),squeeze(d_recon_lanb(:,5,5,:,6))-squeeze(d(:,5,5,:,6)),...
%     squeeze(d_recon_rsvd(:,5,5,:,6))-squeeze(d(:,5,5,:,6))]);


%% from Matlab to Madagascar
rsf_create(clean,size(reshape(d,100,10^4))');
rsf_write(reshape(d,100,10^4),clean);

rsf_create(noisy,size(reshape(dn,100,10^4))');
rsf_write(reshape(dn,100,10^4),noisy);

rsf_create(decimated,size(reshape(dn0,100,10^4))');
rsf_write(reshape(dn0,100,10^4),decimated);

rsf_create(lmafit,size(reshape(d_recon_lmafit,100,10^4))');
rsf_write(reshape(d_recon_lmafit,100,10^4),lmafit);

rsf_create(svd,size(reshape(d_recon_svd,100,10^4))');
rsf_write(reshape(d_recon_svd,100,10^4),svd);

rsf_create(lanb,size(reshape(d_recon_lanb,100,10^4))');
rsf_write(reshape(d_recon_lanb,100,10^4),lanb);

rsf_create(rsvd,size(reshape(d_recon_rsvd,100,10^4))');
rsf_write(reshape(d_recon_rsvd,100,10^4),rsvd);
