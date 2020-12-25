function [dout]=P_Svd(din,N,K)
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


% Rank reduction on the block Hankel matrix


%      [U,D,V]=svds(din,N); % a little bit slower for small matrix
%      dout=U*D*V';
% %      
    [U,D,V]=svd(din);
%     for j=1:N
%         D(j,j)=D(j,j)*(1-D(N+1,N+1)^K/D(j,j)^K);
%     end    
    
    dout=U(:,1:N)*D(1:N,1:N)*(V(:,1:N)');

return
end