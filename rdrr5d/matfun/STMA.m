function [dout]=STMA(Z, K, e, ws)
% Author      : Oboue et al. 2020
%               Zhejiang University
% 
% Date        : January, 2021
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) Oboue et al. 2020

% stma : soft thresholding moving-average operator
% Z    : input signal 
% K    : decrease factor for cooling lambda
% e    : rational tranfer function coefficient 
% ws   : windows size
% dout : approximation signal
%
%  KEY REFERENCE
%  Oboue et al., 2021, Robust damped rank-reduction method for simultaneous denoising and reconstruction of 5-D seismic data, Geophysics, 86, V71–V89.

D=Z;
x=Z(:);

sizeZ = size(Z);
T = randperm(prod(sizeZ));

IDX = T(1:round(0.999*prod(sizeZ)));

M = opRestriction(prod(sizeZ), IDX);

y = M(x,1);

lambda =K*(norm(max(abs(y)))); 

alpha= std(y)*(norm(max(abs(y))) + sqrt(length(y)));

thr=lambda/alpha;

s = SoftTh(D,thr); 

b = (1/ws)*ones(1,ws);

dout = filter(b,e,s);
end
