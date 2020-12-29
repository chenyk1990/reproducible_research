function SYN1_N(noise,nssak1,nsstk1,nssak3,nsstk3,noises)
% Author      : Yangkang Chen
%               Zhejiang University
%         
% Date        : Dec, 2020
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Yangkang Chen
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

%% adding noise
% randn('state',201314);
% var=0.2;
% dn=d+var*randn(size(d));



d=zeros(401,200);
%dn=zeros(512,50);
%rsf_read(d,'../syn1/2Dsynth-clean.rsf');
randn('state',201314);
% d=yc_scale(d,2);
var=1;
%rsf_read(dn1,'../syn1/syn.rsf');
n=yc_bp(var*randn(size(d)),0.004,0,5,30,40);
dn=d+n;


flow=0;fhigh=40;dt=0.004;N=1;verb=0;
tic
n1=fxymssa_grouse(n,flow,fhigh,dt,N,1,verb);
toc
%figure;imagesc([n,n1,n-n1]);
%colormap(seis);
%snrcyk(d,d1,2)
%about 12s
%snr:4.49

tic
n2=fxymssa_grouse(n,flow,fhigh,dt,N,2,verb);
toc
%figure;imagesc([n,n2,n-n2]);
%colormap(seis);
%snrcyk(d,d2,2)
%about 6s
%about 6.46

N=3;
tic
n11=fxymssa_grouse(n,flow,fhigh,dt,N,1,verb);
toc
%figure;imagesc([n,n1,n-n1]);
%colormap(seis);
%snrcyk(d,d1,2)
%about 12s
%snr:4.49

tic
n22=fxymssa_grouse(n,flow,fhigh,dt,N,2,verb);
toc
%figure;imagesc([n,n2,n-n2]);
%colormap(seis);
%snrcyk(d,d2,2)
%about 6s
%about 6.46

e0=sum(sum(abs(n)))
e1=sum(sum(abs(n1)))
e2=sum(sum(abs(n2)))

tic
n3=fxymssa_grouse(n,flow,fhigh,dt,N,3,verb);
toc
%figure;imagesc([n,n3,n-n3]);
%colormap(seis);
%snrcyk(d,d3,2)
%about 1.5s
%SNR: 5.28 (error too much)


tic
n4=yc_fkt(n,'ps',20);
toc
%figure;imagesc([n,n4,n-n4]);
%colormap(seis);
%snrcyk(d,d4,2)
%about 0.05s
%SNR: 4.3




%size:401*200, 4.698,4.605,0.405,0.02

%[3.1446e+04,9.2834e+03,4.2836e+03] k=3
k=[1,2,3,4,5,6,7,8,9,10];
e0=3.1446e+04*ones(1,10);
e1=[5.6612e+03,7.7556e+03,9.2834e+03,1.0569e+04,1.1671e+04,1.2616e+04,1.3521e+04,1.4370e+04,1.5125e+04,1.5848e+04];
e2=[2.3449e+03,3.4170e+03,4.2836e+03,5.0828e+03,5.8136e+03,6.4654e+03,7.0652e+03,7.6368e+03,8.1701e+03,8.6599e+03];

nn=[k',e0',e1',e2'];


rsf_create(noise,size(n)');
rsf_write(n,noise);

rsf_create(nssak1,size(n1)');
rsf_write(n1,nssak1);

rsf_create(nsstk1,size(n2)');
rsf_write(n2,nsstk1);

rsf_create(nssak3,size(n11)');
rsf_write(n11,nssak3);

rsf_create(nsstk3,size(n22)');
rsf_write(n22,nsstk3);

rsf_create(noises,size(nn)');
rsf_write(nn,noises);


