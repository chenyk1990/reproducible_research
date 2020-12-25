function curve(clean,noisy,decimated,lmafit,svd,lanb,rsvd)
%% data preparition and parameterization 
addpath(genpath('../matfun'));
load ../../drr5d/matfun/yc_hyper5d.mat
[nt,nhx,nhy,nx,ny]=size(hyper5d);
d=hyper5d;

dt=0.004;
nf=2^nextpow2(nt);
flow=0;
fhigh=125;
N1=12; % rank for svd and lmafit
N2=12; % rank for lanb
N3=12; % rank for rsvd
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
mask=genmask(reshape(dn,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
dn0=dn.*mask;

%% data processing - 5D MSSA
tic
[d_recon_lmafit] = mssa5d_lmafit(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N1,K,eps,Niter,a,mask);
toc
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_lmafit,nt*nhx*nhy,nx,ny),2)

tic
[d_recon_svd] = mssa5d_svd(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N1,K,eps,Niter,a,mask);
toc
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_svd,nt*nhx*nhy,nx,ny),2)

tic
[d_recon_lanb] = mssa5d_lanb(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N2,K,eps,Niter,a,mask);
toc
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_lanb,nt*nhx*nhy,nx,ny),2)

tic
[d_recon_rsvd] = mssa5d_rsvd(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N3,K,P,eps,Niter,a,mask);
toc
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(dn,nt*nhx*nhy,nx,ny),2)
snrcyk(reshape(d,nt*nhx*nhy,nx,ny),reshape(d_recon_rsvd,nt*nhx*nhy,nx,ny),2)

% figure;
% imagesc([squeeze(d(:,5,5,:,6)),squeeze(dn(:,5,5,:,6)),squeeze(d_recon_lmafit(:,5,5,:,6)),...
%     squeeze(d_recon_lanb(:,5,5,:,6)),squeeze(d_recon_rsvd(:,5,5,:,6)),...
%     squeeze(d_recon_lmafit(:,5,5,:,6))-squeeze(d(:,5,5,:,6)),squeeze(d_recon_lanb(:,5,5,:,6))-squeeze(d(:,5,5,:,6)),...
%     squeeze(d_recon_rsvd(:,5,5,:,6))-squeeze(d(:,5,5,:,6))]);


%% from Matlab to Madagascar
rsf_create(clean,size(reshape(d,101,32*32*5*5))');
rsf_write(reshape(d,101,32*32*5*5),clean);

rsf_create(noisy,size(reshape(dn,101,32*32*5*5))');
rsf_write(reshape(dn,101,32*32*5*5),noisy);

rsf_create(decimated,size(reshape(dn0,101,32*32*5*5))');
rsf_write(reshape(dn0,101,32*32*5*5),decimated);

rsf_create(lmafit,size(reshape(d_recon_lmafit,101,32*32*5*5))');
rsf_write(reshape(d_recon_lmafit,101,32*32*5*5),lmafit);

rsf_create(svd,size(reshape(d_recon_svd,101,32*32*5*5))');
rsf_write(reshape(d_recon_svd,101,32*32*5*5),svd);

rsf_create(lanb,size(reshape(d_recon_lanb,101,32*32*5*5))');
rsf_write(reshape(d_recon_lanb,101,32*32*5*5),lanb);

rsf_create(rsvd,size(reshape(d_recon_rsvd,101,32*32*5*5))');
rsf_write(reshape(d_recon_rsvd,101,32*32*5*5),rsvd);
