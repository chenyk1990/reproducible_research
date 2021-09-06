function Create_deltav(deltav,Nevents,Nz)
% Author      : Yangkang Chen
% Date        : Aug 2021
% 
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2021 Yangkang Chen
%  
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% 

% Nevents=100000;
Nevents=10000;Nz=101;

randn('state',20202122);
dv=randn(Nz,Nevents);
dv=0.5*yc_scale(yc_meanf(dv,10,1,1),1);
% figure;imagesc(dv);colorbar;

rsf_create(deltav,size(dv)');
rsf_write(dv,deltav);






