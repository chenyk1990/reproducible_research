function LDRR(din,d_lrr,d_lrra,n1,n2,n3,dt,lf,hf,N,NN,n1win,n2win,n3win,r1,r2,r3,ifdamp,verb)
% Author      : Yangkang Chen
% Date        : March, 2020

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
d=zeros(n1,n2*n3);
rsf_read(d,din)
d=reshape(d,n1,n2,n3);

%%% Main program goes here !

if ~ifdamp
    %% LRR
    d1=fxymssa_win(d,lf,hf,dt,N,verb,n1win,n2win,n3win,r1,r2,r3);
    %% LDRR
    d2=fxymssa_win_auto(d,lf,hf,dt,N,verb,n1win,n2win,n3win,r1,r2,r3,2);
    
    d1=reshape(d1,n1,n2,n3);
    d2=reshape(d2,n1,n2,n3);
else
    %% LRRA
    d3=fxydmssa_win(d,lf,hf,dt,N,NN,verb,n1win,n2win,n3win,r1,r2,r3);
    %% LDRRA
    d4=fxydmssa_win_auto(d,lf,hf,dt,N,NN,verb,n1win,n2win,n3win,r1,r2,r3,2);

    d3=reshape(d3,n1,n2,n3);
    d4=reshape(d4,n1,n2,n3);
end

%% from Matlab to Madagascar
if ~ifdamp
    rsf_create(d_lrr,size(d1)');
    rsf_write(d1,d_lrr);
    
    rsf_create(d_lrra,size(d2)');
    rsf_write(d2,d_lrra);
else
    rsf_create(d_lrr,size(d3)');
    rsf_write(d3,d_lrr);
    
    rsf_create(d_lrra,size(d4)');
    rsf_write(d4,d_lrra);
end

return

