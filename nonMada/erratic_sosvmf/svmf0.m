function [ D1,win_len ] = svmf0(D,nfw,ifb,axis,l1,l2,l3,l4)
%MFCYK: median filter along first or second axis for 2D profile
%  IN   D:   	intput data 
%       nfw:    window size
%       ifb:    if use padded boundary (if not, zero will be padded)
%       axis:    temporal sampling interval
%      
%  OUT   D1:  	output data
% 
%  Copyright (C) 2019 Zhejiang University
%  Copyright (C) 2019 Yangkang Chen
%
% Example: dsp/test_svmf0.m
%          dsp/test_sosvmf.m
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

if nargin==0
 error('Input data must be provided!');
end
[n1,n2]=size(D);



if nargin==1
 nfw=7;
 ifb=1;
 axis=2;
end;

if nargin==2
 ifb=1;
 axis=2;    
end

Dtmp=yc_mf(D,nfw,ifb,axis);
medianv=sum(abs(Dtmp(:)))/(n1*n2);

if nargin==4
   l1=2;
   l2=0;
   l3=2;
   l4=4;
end

% nfw should be odd
if mod(nfw,2)==0
    nfw=nfw+1;
end


% calculate length
win_len=zeros(size(D));
for i2=1:n2
   for i1=1:n1
       if abs(Dtmp(i1,i2)) <medianv
          if abs(Dtmp(i1,i2))<medianv/2 
             win_len(i1,i2)=nfw+l1;
          else
             win_len(i1,i2)=nfw+l2; 
          end
       else
           if abs(Dtmp(i1,i2))>medianv*2 
             win_len(i1,i2)=nfw-l4;
          else
             win_len(i1,i2)=nfw-l3; 
          end          
       end
       
   end 
end

if axis==2
   D=D.'; 
   win_len=win_len.';
end
win_len2=(win_len-1)/2;

[n1,n2]=size(D);
nfw2=(nfw-1)/2;

nfw_b=(max([nfw+l1,nfw+l2])-1)/2;
if ifb==1
    D=[flipud(D(1:nfw_b,:));D;flipud(D(n1-nfw_b+1:n1,:))];
else
    D=[zeros(nfw_b,n2);D;zeros(nfw_b,n2)];    
end

% output data
D1=zeros(n1,n2);
for i2=1:n2
   for i1=1:n1
      D1(i1,i2)=median(D(i1+nfw_b-win_len2(i1,i2):i1+nfw_b+win_len2(i1,i2),i2)); 
   end 
end

win_len=win_len2*2+1;
if axis==2
    D1=D1.';
    win_len=win_len.';    
end
return
