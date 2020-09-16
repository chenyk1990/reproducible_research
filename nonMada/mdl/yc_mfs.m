function [ dout ] = yc_mfs(din,nfw,ifb,axis,ntimes)
%YCMFS: median filter along first or second axis for 2D profile
%  IN   D:   	intput data 
%       nfw:    window size
%       ifb:    if use padded boundary (if not, zero will be padded)
%       axis:    temporal sampling interval
%       tmes:   number of filtering times
%      
%  OUT   D1:  	output data
% 
%  Copyright (C) 2018 The University of Texas at Austin
%  Copyright (C) 2018 Yangkang Chen
%
% Example: dsp/test_yc_meanf_repeat.m
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
% REFERENCE
% For more details about median filtering in seismics, please refer to 
% Chen, Y., J. Yuan, Z. Jin, K. Chen and L. Zhang, 2014, Deblending using normal moveout and median filtering in common-midpoint gather, Journal of Geophysics and Engineering, 11, 045012.
% Chen, Y., 2015, Deblending using a space-varying median filter, Exploration Geophysics, 46, 332-341.
% Gan, S., S. Wang, Y. Chen, X. Chen, and Kui Xiang, 2016, Separation of simultaneous sources using a structural-oriented median filter in the flattened dimension, Computers & Geosciences, 86, 46-54.
% Huang, W., R. Wang, X. Gong, and Y. Chen, 2018, Iterative deblending of simultaneous- source seismic data with structuring median constraint, IEEE Geoscience and Remote Sensing Letters, 15, 58-62.
% Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805?1823.
% Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2020, Erratic-noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting


if nargin==0
 error('Input data must be provided!');
end

if nargin==1
 nfw=7;
 ifb=1;
 axis=2;
end

if nargin==2
 ifb=1;
 axis=2;    
end

dtmp=din;
if axis==3
for itime=1:ntimes
dtmp=yc_mf(dtmp,nfw,ifb,1);
dtmp=yc_mf(dtmp,nfw,ifb,2);
end    
else
for itime=1:ntimes
dtmp=yc_mf(dtmp,nfw,ifb,axis);
end
end
dout=dtmp;


return
