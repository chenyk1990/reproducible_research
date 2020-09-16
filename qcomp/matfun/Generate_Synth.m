 function Generate_Synth(d_c,d_a,d_n)
% Author      : Wang et al.
%               Zhejiang Universit
%         
% Date        : Jan, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Wang et al. 
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

%%from Madagascar to Matlab
% create memory

n1=1001;
n2=12;
dt=0.002;
t=0:dt:(n1-1)*dt;
x=1:n2;
df=1/(n1*dt);
freqdom=30;
Q = 50*ones(1,n1);

%% wavelet and wavelet matrix

Wavn1=n1;
j=sqrt(-1);
tw=-(Wavn1-1)/2*dt:dt:(Wavn1-1)/2*dt;
wavelet=(1-2*(pi*freqdom.*tw).^2).*exp(-(pi*freqdom.*tw).^2);
wavelet=wavelet./max(abs(wavelet));
% figure;
% plot(wavelet);

MatrixWavelet=zeros(n1,n1);
nn=ceil(n1/2);
for i=1:nn
    MatrixWavelet(1:nn+i-1,i)=wavelet(nn-i+1:end);
end
for i=1:nn-1
    MatrixWavelet(i+1:end,i+nn)=wavelet(1:(end-i));
end
% figure;
% yc_wigb(MatrixWavelet);


%% reflectivity
reflect=zeros(n1,n2);

for i=1:n2
    reflect(100,i)=0.03;
    reflect(160+6*i,i)=-0.03;
    reflect(260-8*i,i)=0.03;
    reflect(390,i)=0.03;
    reflect(400+4*i,i)=-0.03;
    reflect(600,i)=0.03;
    reflect(700+10*i,i)=-0.03;
end

%% signal
signal=MatrixWavelet*reflect;
Signal=zeros(n1,n2);
for i=1:n2
    for k=1:n1
        Signal(k,i)=signal(k+(i-1)*n1);
    end
end

%% time-domain absorption matrix
[MatrixAbosrp_fre,MatrixAbosrp]=time_absorp_wavelet_matrix(n1,dt,freqdom,Q);

%% attenuated signal
Signal_absorp=MatrixAbosrp*Signal;

fp=fopen('Signal_absorp_noisy.dat','rb');
Signal_absorp_noisy=fread(fp,[n1, n2],'float');
fclose(fp);

dc=Signal;
da=Signal_absorp;
dn=Signal_absorp_noisy;

yc_snr(dc,dn)

rsf_create(d_c,size(dc)');
rsf_write(dc,d_c);

rsf_create(d_a,size(da)');
rsf_write(da,d_a);

rsf_create(d_n,size(dn)');
rsf_write(dn,d_n);

return
