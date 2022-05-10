function [ m,d ] = yc_adjnull( adj,add,nm,nd,m,d )
%% Claerbout-style adjoint zeroing Zeros out the output (unless add is true). 
% Useful first step for and linear operator.
%  
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) and later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%%
% adj : adjoint flag; add: addition flag; nm: size of m; nd: size of d
if(add)
    return
end

if(adj)
    m=zeros(nm,1);
    for i = 0 : nm-1
        m(i+1) = 0.0;
    end
else
    d=zeros(nd,1);
    for i = 0 : nd-1
        d(i+1) = 0.0;
    end    

end




