function [ simi ] = localsimi(d1,d2,rect,niter,eps,verb)
%  LOCALSIMI: calculate local similarity between two datasets
%
%  IN   d1:   	input data 1
%       d2:     input data 2
%       verb:   verbosity flag (default: 0)
%
%  OUT  simi:  	calculated local similarity, which is of the same size as d1 and d2
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
%  Reference:   
%  				1. Chen, Y. and S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, , 80, WD1-WD9. (Note that when local similarity is used for noise suppression purposes, this reference must be cited.)
%               2. Local seismic attributes, Fomel, Geophysics, 2007
%
% DEMO
% test_localortho.m 

if nargin<2
    error('Input data 1 and data 2 must be provided!');
end

[n1,n2,n3]=size(d1);

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
end;

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
[ ratio ] = str_divne(d2, d1, niter, rect, ndat, eps_dv, eps_cg, tol_cg,verb);
[ ratio1 ] = str_divne(d1, d2, niter, rect, ndat, eps_dv, eps_cg, tol_cg,verb);

simi=sqrt(abs(ratio.*ratio1));






