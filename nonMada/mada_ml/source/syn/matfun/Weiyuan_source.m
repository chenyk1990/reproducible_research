function Weiyuan_source(slocname,Nevents)
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

xm=[0,4.0];
ym=[0,3.3];
zm=[0,3.8];

%% purely random (not on the grid)
rand('state',202122);
x=rand(Nevents,1)*(xm(2)-xm(1))+xm(1);%scale transformation
rand('state',20212223);
y=rand(Nevents,1)*(ym(2)-ym(1))+ym(1);%scale transformation
rand('state',2021222324);
z=rand(Nevents,1)*(zm(2)-zm(1))+zm(1);%scale transformation

% %% random but on the grid
% x=[0:0.05:5]';nx=length(x);
% y=[0:0.05:5]';ny=length(y);
% z=[0:0.05:5]';nz=length(z);

x=[0:0.0331:100*0.0331]';nx=length(x);
y=[0:0.03835:100*0.03835]';ny=length(y);
z=[0:0.04:100*0.04]';nz=length(z);

rand('state',202122);
x=x(floor(rand(Nevents,1)*nx)+1);
rand('state',20212223);
y=y(floor(rand(Nevents,1)*ny)+1);
rand('state',20212224);
z=z(floor(rand(Nevents,1)*nz)+1);

% x=x*0+3.5;
% y=y*0+2.5;
% z=z*0+2.5;
%% create locations
fid=fopen(slocname,'w');
fprintf(fid,'%3.4f %3.4f %3.4f\n',[x,y,z]');

%figure;
%scatter3(x,y,z);



