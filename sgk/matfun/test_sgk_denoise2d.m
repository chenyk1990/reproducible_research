clc;clear;close all;
%% DEMO script for DL (SGK) based denoising (2D)
%
% Key reference:
% Chen, Y., 2017, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 209, 21-31.
% Chen, Y., W. Huang, D. Zhang, and W. Chen, 2016, An open-source Matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% 
% More References:
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2019, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 218, 1379?1397. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51?KS61.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.
% etc. 

ne=20;%number of events;
dt = 4./1000;
tmax = dt*127;
h = [-160:10:150];
tau=linspace(0.1,1.8,ne);
v0=linspace(1800,3000,ne);%exact
randn('state',201819202122);amp = randn(ne,1);%[1., -1.,1];
f0 = 20;
snr = 2000;%default snr=2
L = 6;
seed=201517;
dc=hevents(dt,f0,tmax,h,tau,v0,amp,snr,L,seed);dc=yc_scale(dc,2);

[n1,n2]=size(dc);dt=0.004;
t=[0:n1-1]*dt;x=1:n2;
%%
randn('state',201617);
d=dc+0.1*randn(size(dc));
[n1,n2]=size(d);

%% patch size l1*l2
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
param.T=1;      %sparsity level
% param.D=DCT;    %initial D
param.niter=10; %number of SGK iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
%param.exact:   Exact SGK update or approximate
param.K=64;     %number of atoms, dictionary size
%for X=DG
%size of X: MxN
%size of D: MxK
%size of G: KxN

%% Option 1: denoise only using the integrated function
param=struct('T',2,'niter',10,'mode',1,'K',64,'D',DCT);param=rmfield(param,'D');
mode=1;l1=4;l2=4;s1=2;s2=2;perc=100;
tic
[d1]=yc_sgk_denoise(d,mode,[l1,l2,1],[s1,s2,1],perc,param);
toc
figure;imagesc([dc,d,d1,d-d1]);colormap(jet);
yc_snr(dc,d1)
%11.8746 dB, 0.4s
%% benchmark with KSVD (worse but much faster)
tic
d2=yc_ksvd_denoise(d,mode,[l1,l2,1],[s1,s2,1],perc,param);
toc
figure;imagesc([dc,d,d2,d-d2]);colormap(jet);
yc_snr(dc,d2)
%11.90 dB, 2.8s

%% Option 2: without initialized param.D
% param=rmfield(param,'D');
% param=struct('T',3,'niter',10,'mode',1,'K',64);
% mode=1;l1=8;l2=8;s1=4;s2=4;perc=7;
% d1=yc_sgk_denoise(d,mode,[l1,l2,1],[s1,s2,1],perc,param);
% figure;imagesc([dc,d,d1,d-d1]);colormap(seis);
% yc_snr(dc,d1)



