function BIN(datain,x0,y0,dataout,nt,ntraces)
% Author      : Yangkang Chen
% Date        : Mar 2021
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
%nt=1001
%ntraces=918

x00=zeros(ntraces,1);
y00=zeros(ntraces,1);
data=zeros(nt,ntraces);
rsf_read(x00,x0);
rsf_read(y00,y0);
rsf_read(data,datain); 

x0=x00*0.3048;
y0=y00*0.3048;

%% binning (put a 2D data into a 3D volume)
% about 100 receivers per arm
nx=201;
ny=201;
ox=0;
oy=0;
mx=ox+200*24;
my=oy+200*24;
dt=0.002;
[ dout,xout,yout,mask ] = yc_bin3d(data,x0,y0,nx,ny,ox,oy,mx,my);


%% from Matlab to Madagascar
dout=reshape(dout,nt,nx*ny);
rsf_create(dataout,size(dout)');
rsf_write(dout,dataout);







