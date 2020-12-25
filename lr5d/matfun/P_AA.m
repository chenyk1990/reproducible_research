function [dout]=P_AA(din,nf,nhx,nhy,nx,ny,lx,ly,lxx,lhx,lhy,lhxx,lhyy)
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



% Averaging the block Hankel matrix to output the result
dout=zeros(nhx,nhy,nx,ny);
M2=P_A2(din,nhx,nhy,lx,ly,lxx,lhx,lhy,lhxx,lhyy);

for i=1:nhx
    for j=1:nhy
        dm=reshape(M2(i,j,:,:),lhx*lhy,lhxx*lhyy);
        dori=P_A(dm,nx,ny,lhx,lhy);
        dout(i,j,:,:)=dori;
    end
end
return
end

function [dout]=P_A2(din,nhx,nhy,lx,ly,lxx,lhx,lhy,lhxx,lhyy)
% Averaging the block Hankel matrix to output the result
dout=zeros(nhx,nhy,lhx*lhy,lhxx*lhyy);

for j=1:nhy
    if j<ly
        for id=1:j
            dout(:,j,:,:) =dout(:,j,:,:)+ ave_antid2(din(1+(j-1)*lx*lhx*lhy-...
                (id-1)*lx*lhx*lhy:j*lx*lhx*lhy-(id-1)*lx*lhx*lhy,1+(id-1)*lxx*lhxx*lhyy:lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy),nhx,lx,lhx,lhy,lhxx,lhyy)/j;
        end
    else
        for id=1:(nhy-j+1)           
            dout(:,j,:,:) =dout(:,j,:,:)+ ave_antid2(din((ly-1)*lx*lhx*lhy+...
                1-(id-1)*lx*lhx*lhy:ly*lx*lhx*lhy-(id-1)*lx*lhx*lhy,(j-ly)*lxx*lhxx*lhyy+1+(id-1)*lxx*lhxx*lhyy:(j-ly+1)*lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy),nhx,lx,lhx,lhy,lhxx,lhyy)/(nhy-j+1);
        end
    end
end
return
end

function [dout]=P_A(din,nhx,nhy,lx,ly)
% Averaging the block Hankel matrix to output the result
lxx=nhx-lx+1;
lyy=nhy-ly+1;
dout=zeros(nhx,nhy);

for j=1:nhy
    if j<ly
        for id=1:j
            dout(:,j) =dout(:,j)+ ave_antid(din(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx))/j;
        end
    else
        for id=1:(nhy-j+1)
            dout(:,j) =dout(:,j)+ ave_antid(din((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx))/(nhy-j+1);
        end
    end
end
return
end

function [dout] =ave_antid2(din,nhx,lx,lhx,lhy,lhxx,lhyy);
% averaging along antidiagonals
dout=zeros(nhx,lhx*lhy,lhxx*lhyy);

for i=1:nhx
    if i<lx
        for id=1:i            
            dout(i,:,:)=dout(i,:,:) + reshape(din(1+(i-id)*lhx*lhy:(i-id+1)*lhx*lhy,1+(id-1)*lhxx*lhyy:id*lhxx*lhyy)/i,1,lhx*lhy,lhxx*lhyy);
        end
    else
        for id=1:nhx+1-i
            dout(i,:,:)=dout(i,:,:) + reshape(din(1+(lx-id)*lhx*lhy:(lx-id+1)*lhx*lhy,1+(i-lx+id-1)*lhxx*lhyy:(i-lx+id)*lhxx*lhyy)/(nhx+1-i),1,lhx*lhy,lhxx*lhyy);
        end
    end
end
dout=reshape(dout,nhx,1,lhx*lhy,lhxx*lhyy);
return
end

function [dout] =ave_antid(din);
% averaging along antidiagonals
[n1,n2]=size(din);
nout=n1+n2-1;
dout=zeros(nout,1);
for i=1:nout
    if i<n1
        for id=1:i
            dout(i)=dout(i) + din(i-(id-1),1+(id-1))/i;
        end
    else
        for id=1:nout+1-i
            dout(i)=dout(i) + din(n1-(id-1),1+(i-n1)+(id-1))/(nout+1-i);
        end
    end
end
return
end

