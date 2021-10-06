function Real(noisy,denoised1,denoised2,n1,n2,n3,T,niter,K,ll1,ll2,ll3,ss1,ss2,ss3,perc)
% Author      : Hang Wang and Yangkang Chen
%               Zhejiang University
%         
% Date        : April, 2020
%
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%
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
%  Reference:
%  Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE TGRS, doi: 10.1109/TGRS.2020.3030740
%  Chen, 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, GJI, 222, 1717-1727.

dn=zeros(n1,n2*n3);
d1=zeros(n1,n2*n3);

rsf_read(dn,noisy);dn=reshape(dn,n1,n2,n3);
rsf_read(d1,denoised1);d1=reshape(d1,n1,n2,n3);


%% simultaneous denoising and reconstruction
%% SGK
param=struct('T',T,'niter',niter,'mode',1,'K',K);
mode=1;
l1=ll1;
l2=ll2;
l3=ll3;
s1=ss1;
s2=ss2;
s3=ss3;
perc=perc;

%% create initial DCT
c1=l1;c2=l2;c3=l3;
Dinit=yc_initD([l1,l2,l3],[c1,c2,c3]);
Dinit=Dinit(:,1:K);
param.D=Dinit;

% %%test yc_initD
% l1=4;l2=4;l3=1;
% c1=l1;c2=l2;c3=l3;
% Dinit=yc_initD([l1,l2,l3],[c1,c2,c3]);
% % figure;imagesc(Dinit);

if n3==1
XX=yc_patch(d1,mode,l1,l2,s1,s2);
XXn=yc_patch(dn-d1,mode,l1,l2,s1,s2);
else
XX=yc_patch3d(d1,mode,l1,l2,l3,s1,s2,s3);
XXn=yc_patch3d(dn-d1,mode,l1,l2,l3,s1,s2,s3);   
end
tic
[DD,GG]=yc_sgk(XX,param);
toc
tic
[DDksvd,GGksvd]=yc_ksvd(XX,param);
toc
Gn=yc_ompN(DD,XXn,T);
Gn=yc_pthresh(Gn,'ph',perc);
Xn=DD*Gn;
if n3==1
d11=yc_patch_inv(Xn,mode,n1,n2,l1,l2,s1,s2);
else
d11=yc_patch3d_inv(Xn,mode,n1,n2,n3,l1,l2,l3,s1,s2,s3);
end
d2=d1+d11;

d2=reshape(d2,n1,n2*n3);
rsf_create(denoised2,size(d2)');
rsf_write(d2,denoised2);

























