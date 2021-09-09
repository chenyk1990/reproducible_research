function [ signal2,noise2,low ] = localortho(signal,noise,rect,niter,eps,verb)
%  YC_LOCAORTHO: Noise attenuation using local signal-and-noise
%  orthogonalization and output the local orthogonalization weight (LOW)
%
%  IN   signal:    initial signal
%       noise:     initial noise
%       rect:   3-D vector denoting smooth radius
%       niter:  number of CG iterations
%       eps:    regularization parameter, default 0.0
%       verb:   verbosity flag (default: 0)
%
%  OUT  signal2: orthogonalized signal
%       noise2:  orthogonalized noise
%       low:   local orthogonalization weight
%
%  Copyright (C) 2016 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  Reference:   1. Random noise attenuation using local signal-and-noise orthogonalization
%               Chen and Fomel, 2015, Geophysics
%               2. Ground-Roll Noise Attenuation Using a Simple and Effective Approach Based on 
%               Local Band-Limited Orthogonalization, Chen et al., 2015, IEEE Geoscience and Remote Sensing Letters
%               3. Iterative deblending with multiple constraints based on shaping regularization,
%               Chen, 2015, IEEE Geoscience and Remote Sensing Letters
%               4. Orthogonalized morphological reconstruction for weak signal detection in micro-seismic monitoring:
%               Methodology, Huang et al., 2018, GJI
%               5. Surface-related multiple leakage extraction using local primary-and-multiple 
%               orthogonalization, Zhang et al., 2020, Geophysics
%               6. Non-stationary local signal-and-noise orthogonalization, Chen et al.,
%               2020, Geophysics
%               7. Local primary-and-multiple orthogonalization for leaked internal multiple crosstalk estimation and attenuation on full-wavefield migrated images
%               Zhang, et al., 2020, Geophysics
%  Example:
%  test_localortho.m
% 
if nargin<2
    error('Input data 1 and data 2 must be provided!');
end

[n1,n2,n3]=size(signal);

if nargin==2
    rect=ones(3,1);
    if(n1==1) error('data must be a vector or a matrix!');
    else
        rect(1)=20;
    end
    if(n2~=1) rect(2)=10;end
    if(n3~=1) rect(3)=10;end
    niter=50;
    eps=0.0;    
    verb=1;
end

if nargin==3
   niter=50;
   eps=0.0;
   verb=1;
end

if nargin==4
   eps=0.0;
   verb=1;
end

if nargin==5
   verb=1;
end

%eps=0.0;
nd=n1*n2*n3;
ndat=[n1,n2,n3];

eps_dv=eps;
eps_cg=0.1; 
tol_cg=0.000001;
[ ratio ] = str_divne(noise, signal, niter, rect, ndat, eps_dv, eps_cg, tol_cg,verb);

signal2=signal+ratio.*signal;
noise2=noise-ratio.*signal;
low=ratio;











