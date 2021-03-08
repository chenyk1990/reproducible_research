function [dout]=atanTma1(y, ws, e1, o1, b1)
% atanTma1 : Improved proximity function by mixing the moving-average filter and the arctangent penalty function.
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
%       e1:     rational tranfer function coefficient
%       o1:     decrease factor for cooling lambda
%       b1:     penalty convexity parameter 
% OUTPUT : approximation signal

%  KEY REFERENCE
%  Oboue and Chen 2021, Enhanced low-rank matrix estimation for simultaneous denoising and reconstruction of five-dimensional seismic data
% 

c = (1/ws)*ones(1,ws);

x = filter(c,e1,y);

lambda =o1*(norm(max(abs(y(:))))); 

alpha= std(y(:))*(norm(max(abs(y(:)))) + sqrt(length(y(:))));

T=lambda/alpha; % threshold (scalar or multidimensional array)

dout = atanT(x, T, b1);
end
