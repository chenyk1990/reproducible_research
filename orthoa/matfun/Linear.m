 function Linear(fin,fxp,fx,n1,n2,N,mu,verb)
% Author      : Guangtan Huang and Yangkang Chen
%         
% Date        : Oct, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Guangtan Huang and Yangkang Chen
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
din=zeros(n1,n2);
rsf_read(din,fin);

%% Main program goes here !
flow=0;
fhigh=120;
dt=0.004;
lf=N;
mu=mu;
verb=verb;
d1 = fxarma(din,flow,fhigh,dt,lf,mu,verb);
d2 = fxdecon(din,flow,fhigh,dt,lf,mu,verb);

%% from Matlab to Madagascar
rsf_create(fxp,size(d1)');
rsf_write(d1,fxp);

rsf_create(fx,size(d2)');
rsf_write(d2,fx);

