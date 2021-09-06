function [ dout,xout,yout,mask ] = yc_bin3d(din,x,y,nx,ny,ox,oy,mx,my)
%yc_bin3d (previous cykbin3d): 3D seismic data binning (including 1D row vector and 2D seismics)
%  IN    d:   	intput 2D data
%        x:     input x coordinates
%        y:     input y coordinates
%        nx:    input number of binned x points
%        ny:    input number of binned y points
%        ox:    min of x
%        oy:    min of y
%        mx:    max of x
%        my:    max of y
%
%  OUT   d1:  	output data
%        xout:  output x coordinates
%        yout:  output x coordinates
%        mask:  mask operator for interpolation
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
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
%  Demo:
%  ~/chenyk/friends/chengping/matfun/test_inter_rf.m
%  ~/chenyk/friends/chengping/matfun/test_inter_chengping.m
%
%  see also: cykbin2d
%
[n1,n2]=size(din);

if nargin==5
    ox=min(x);
    oy=min(y);
    mx=max(x);
    my=max(y);
    dx=(mx-ox)/(nx-1);
    dy=(my-oy)/(ny-1);
else
    dx=(mx-ox)/(nx-1);
    dy=(my-oy)/(ny-1);
end
xout=ox:dx:mx;
yout=oy:dy:my;

dout=zeros(n1,nx,ny);
mask=ones(nx,ny);

% dx
% max(x)
% min(x)

for iy=1:ny
    for ix=1:nx
        index=find(x>=xout(ix) & x<xout(ix)+dx & y>=yout(iy) & y<yout(iy)+dy);
        n=length(index);
        if n==1;
            dout(:,ix,iy)=din(:,index);
        else if n==0;
                mask(ix,iy)=0;
                dout(:,ix,iy)=zeros(n1,1);
            end
            
            if n>=2
                if x(index(1))==xout(ix) & y(index(1))==yout(iy);
                    dout(:,ix,iy)=din(:,index(1));
                else
                    t1=sqrt((x(index(1))-xout(ix))^2+(y(index(1))-yout(iy))^2);
                    t2=sqrt((x(index(2))-xout(ix))^2+(y(index(2))-yout(iy))^2);
                    dout(:,ix,iy)=(t1*din(:,index(2))+t2*din(:,index(1)))/(t1+t2);
%                     dout(:,ix)=(t1*din(:,index(2))+t2*din(:,index(1)))/(t1+t2);
                end
            end
        end
    end
end

if nargout==4
    mask=ones(n1,1)*reshape(mask,1,nx*ny);
    mask=reshape(mask,n1,nx,ny);
end

return