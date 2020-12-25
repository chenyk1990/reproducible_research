function [dout]=P_HH(din,lx,ly,lxx,lhx,lhy,lhxx,lhyy)
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


H2=P_H2(din,lhx,lhy,lhxx,lhyy);
dout=P_H4(H2,lx,ly,lxx,lhx,lhy,lhxx,lhyy);

return
end

function [dout]=P_H4(din,lx,ly,lxx,lhx,lhy,lhxx,lhyy)
% forming block Hankel matrix
[~,nhy,~,~]=size(din);

for j=1:nhy
    h2=squeeze(din(:,j,:,:));
    r=P_H3(h2,lhx,lhy,lhxx,lhyy);
    if j<ly
        for id=1:j
            dout(1+(j-1)*lx*lhx*lhy-(id-1)*lx*lhx*lhy:j*lx*lhx*lhy-(id-1)*lx*lhx*lhy,1+(id-1)*lxx*lhxx*lhyy:lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy) = r;
        end
    else
        for id=1:(nhy-j+1)
            dout((ly-1)*lx*lhx*lhy+1-(id-1)*lx*lhx*lhy:ly*lx*lhx*lhy-(id-1)*lx*lhx*lhy,(j-ly)*lxx*lhxx*lhyy+1+(id-1)*lxx*lhxx*lhyy:(j-ly+1)*lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy)=r;
        end
    end
end
return
end

function [dout]=P_H3(din,lhx,lhy,lhxx,lhyy)
% forming 3rd order block Hankel matrix

[nhx,~,~]=size(din);
lx=floor(nhx/2)+1;

for j=1:nhx
    if j<lx
        for id=1:j
            dout(1+(j-1)*lhx*lhy-(id-1)*lhx*lhy:j*lhx*lhy-(id-1)*lhx*lhy,1+...
                (id-1)*lhxx*lhyy:lhxx*lhyy+(id-1)*lhxx*lhyy) = din(j,:,:);
        end
    else
        for id=1:(nhx-j+1)
            dout((lx-1)*lhx*lhy+1-(id-1)*lhx*lhy:lx*lhx*lhy-(id-1)*lhx*lhy,...
                (j-lx)*lhxx*lhyy+1+(id-1)*lhxx*lhyy:(j-lx+1)*lhxx*lhyy+(id-1)*lhxx*lhyy)=din(j,:,:);
        end
    end
end

return
end

function [dout]=P_H2(din,lhx,lhy,lhxx,lhyy)
% forming block Hankel matrix
[nhx,nhy,~,~]=size(din);

dout=zeros(nhx,nhy,lhx*lhy,lhxx*lhyy);

for ii=1:nhx
    for jj=1:nhy
        d=squeeze(din(ii,jj,:,:));
        dout(ii,jj,:,:)=P_H(d,lhx,lhy);
    end
end
return
end

function [dout]=P_H(din,lx,ly)
% forming block Hankel matrix
[nhx,nhy]=size(din);
lxx=nhx-lx+1;
lyy=nhy-ly+1;

for j=1:nhy
    r=hankel(din(1:lx,j),[din(lx:nhx,j)]);
    if j<ly
        for id=1:j
            dout(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx) = r;
        end
    else
        for id=1:(nhy-j+1)
            dout((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx)=r;
        end
    end
end
return
end


