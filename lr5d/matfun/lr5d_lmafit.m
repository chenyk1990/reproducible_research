function [d_recon] = lr5d_lmafit(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N,K,eps,Niter,a,mask)
%% LR5D
% 
%  Copyright (C) 2020 Yangkang Chen (with contributions from all authors of Wu et al., 2020)
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details:
%  http://www.gnu.org/licenses/
%
%  References:   
%
%  Wu et al., 2020, Fast and robust low-rank approximation for high-dimensional seismic data reconstruction, IEEE Access, 8, 175501-175512.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.



% Transform into F-X domain
DATA_FX=fft(dn0,nf,1);
DATA_FX0=zeros(nf,nhx,nhy,nx,ny);

MASK=squeeze(mask(1,:,:,:,:));

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1
    ihigh=floor(nf/2)+1;
end

lhx=floor(nx/2)+1;
lhxx=nx-lhx+1;
lhy=floor(ny/2)+1;
lhyy=ny-lhy+1;

lx=floor(nhx/2)+1;
lxx=nhx-lx+1;
ly=floor(nhy/2)+1;
lyy=nhy-ly+1;

% main loop on frequency
for k=ilow:ihigh

    S_obs=squeeze(DATA_FX(k,:,:,:,:)); 
    Sn_1=S_obs;
    
    for iter=1:Niter
        
        M4=P_HH(Sn_1,lx,ly,lxx,lhx,lhy,lhxx,lhyy);
        M4=P_Lmafit(M4,N);
        Sn=P_AA(M4,nf,nhx,nhy,nx,ny,lx,ly,lxx,lhx,lhy,lhxx,lhyy);
        
        Sn=a(iter)*S_obs+(1-a(iter))*MASK.*Sn+(1-MASK).*Sn;
        
        if norm(reshape(Sn,nhx*nhy,nx*ny)-reshape(Sn_1,nhx*nhy,nx*ny),'fro')<eps
            break;
        end
        
        Sn_1=Sn;
        
    end
    
    DATA_FX0(k,:,:,:,:)=reshape(Sn,1,nhx,nhy,nx,ny);
   
    if(mod(k,5)==0)
        fprintf( 'F %d is done!\n\n',k);
    end
    
end

% Back to TX (the output)
for k=nf/2+2:nf
    DATA_FX0(k,:,:,:,:) = conj(DATA_FX0(nf-k+2,:,:,:,:));
end

d_recon=real(ifft(DATA_FX0,[],1));
d_recon=d_recon(1:nt,:,:,:,:);

end