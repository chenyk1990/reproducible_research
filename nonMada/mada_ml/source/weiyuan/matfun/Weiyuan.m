function Weiyuan(dtimes,data,label,location_name,Nevents,Nstations)
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
times=zeros(101,101*Nevents);
rsf_read(times,dtimes);
times=reshape(times,101,101,Nevents);

%station coordinates in [lon,lat];
weiyuan_stas=load('stations_weiyuan.txt');
%station coordinates in [x,y];
weiyuan_xys=[weiyuan_stas(:,1)-min(weiyuan_stas(:,1)),weiyuan_stas(:,2)-min(weiyuan_stas(:,2))]*111.1949;

%axes
x=[0:0.0331:100*0.0331]';nx=length(x);
y=[0:0.03835:100*0.03835]';ny=length(y);
z=[0:0.04:100*0.04]';nz=length(z);

data222=zeros(3,Nstations,Nevents);%[x,y,t] x Nstation x Nevent
% bin onto regular grids
nsta=size(weiyuan_xys,1);
nx=101;
ny=101;
nz=101;
ox=min(weiyuan_xys(:,1));
oy=min(weiyuan_xys(:,2));
mx=max(weiyuan_xys(:,1));
my=max(weiyuan_xys(:,2));
[ dout,xout,yout,mask ] = yc_bin3d(zeros(10,nsta),weiyuan_xys(:,1),weiyuan_xys(:,2),nx,ny,ox,oy,mx,my);
inds_wy=find(mask(1,:,:)==1);
%% only save real station traveltimes
for ie=1:Nevents
    t=times(:,:,ie);
    [tt,inds]=sort(t(inds_wy),'ascend');
    
    for is=1:Nstations
        ix=mod(inds_wy(inds(is)),101);
        if ix==0
            nx=101;ix=nx;
        end
        iy=ceil(inds_wy(inds(is))/101);
    data222(1,is,ie)=x(ix);
    data222(2,is,ie)=y(iy);
    data222(3,is,ie)=tt(is);
    end
end

%location_name='locations_100000.txt';
fid=fopen(location_name,'r');
label_location=fscanf(fid,'%f %f %f\n',[3,Nevents]);


%save ml_source_train_100000_weiyuan4.mat data222 label_location


%% from Matlab to Madagascar
rsf_create(data,size(data222(:,:))');
rsf_write(data222(:,:),data);

rsf_create(label,size(label_location)');
rsf_write(label_location,label);


