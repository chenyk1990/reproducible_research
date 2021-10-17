function Real2(din,dmask,dout,n1,n2,n3,dt,lf,hf,N,Niter,mode)
% Author      : Yangkang Chen
% Date        : Dec, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API

%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Yangkang Chen
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

%%from Madagascar to Matlab
% create memory

d=zeros(n1,n2*n3);
dm=zeros(n1,n2*n3);
rsf_read(d,din)
rsf_read(dm,dmask)
d=reshape(d,n1,n2,n3);
dm=reshape(dm,n1,n2,n3);

%%% Main program goes here !
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
NN=3;
verb=1;
d1 = fxydmssa_recon(d,dm,lf,hf,dt,N,NN,Niter,eps,verb,mode,a);

%% from Matlab to Madagascar
rsf_create(dout,size(d1)');
rsf_write(d1,dout);
    



return

