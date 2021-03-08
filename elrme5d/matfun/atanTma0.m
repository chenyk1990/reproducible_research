function [dout]=atanTma0(y, ws, e0, o0, b0)
% atanTma0 : Improved proximity function by mixing the moving-average filter and the arctangent penalty function.
% Author      : Oboue and Chen, 2021
%               Zhejiang University
% 
% Date        : February, 2021
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2021 Zhejiang University
%  Copyright (C) Oboue and Chen, 2021

% INPUT
%       y : data (scalar or multidimensional array)
%       ws:     windows size 
%       e0:     rational tranfer function coefficient
%       o0:     decrease factor for cooling lambda
%       b0:     penalty convexity parameter 
% OUTPUT : approximation signal

%  KEY REFERENCE
%  Oboue and Chen 2021, Enhanced low-rank matrix estimation for simultaneous denoising and reconstruction of five-dimensional seismic data
% 

c = (1/ws)*ones(1,ws);

x = filter(c,e0,y);

lambda = o0*sqrt(length(x)); 
alpha= length(x);

T=lambda/alpha; % threshold (scalar or multidimensional array)

dout = atanT(x, T, b0);
end
