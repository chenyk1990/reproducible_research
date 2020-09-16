 function Qcdn(d_n,d_d,n1,n2,n3,dt,freqdom,liter,niter,Q0,verb,ifthr,perc)
% Author      : Hang Wang and Yangkang Chen
%         
% Date        : Mar, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Hang Wang and Yangkang Chen
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
% REFERENCE
% Wang et al., 2020, TGRS, Q-compensated denoising of seismic data

%%from Madagascar to Matlab
% create memory
if nargin==10
   ifthr=0;%if secondary curvelet thresholding
   perc=100;
end
if nargin==11
   ifthr=0;%if secondary curvelet thresholding
   perc=100;
end


dn=zeros(n1,n2*n3);
rsf_read(dn,d_n)

t=0:dt:(n1-1)*dt;

%% wavelet and wavelet matrix
j=sqrt(-1);
tw=-(n1-1)/2*dt:dt:(n1-1)/2*dt;
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

Q = Q0*ones(1,n1);
[MatrixAbosrp_fre,MatrixAbosrp]=time_absorp_wavelet_matrix(n1,dt,freqdom,Q);

A = MatrixWavelet*MatrixAbosrp;
Param.W=A;

for tr=1:n2*n3

    b = dn(:,tr);
    
    [xDCA_1,misfit] = yc_pcg(@mtx_oper,Param,b,zeros(n1,1),liter,niter,~verb);
    b_comp_DCA_1 = MatrixWavelet*xDCA_1;
    comp_DCA_1(:,tr)=b_comp_DCA_1;
    r(:,tr)=xDCA_1;
    
    fprintf("i3=%d/%d is done\n",round(tr/n2),n3);
end

d1=comp_DCA_1;

% remove nan
tmp=find(isnan(d1));
d1(tmp)=zeros(size(tmp));

if ifthr
if n3==1
d1=yc_fdct_wrapping_thr(d1,perc); %remove remaining noise
else
d1=reshape(d1,n1,n2,n3);
d1=yc_fdct3d_thr(d1,perc);		  %remove remaining noise
d1=reshape(d1,n1,n2*n3);
end
end

rsf_create(d_d,size(d1)');
rsf_write(d1,d_d);

return
