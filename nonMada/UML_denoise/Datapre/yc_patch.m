function [ X ] = yc_patch( A,mode,l1,l2,o1,o2 )
%decompose the image into patches:
%  
% by Yangkang Chen
% Oct, 2017
%
% Input: 
%   D: input image
%   mode: patching mode
%   l1: first patch size
%   l2: second patch size
%   o1: first shifting size
%   o2: second shifting size
%   
% Output:
%   X: patches
% 
% Modified on Dec 12, 2018 (the edge issue, arbitrary size for the matrix)
% 
% Examples:
%    ~/chenyk/published/sgk/matfun/demo_omp.m

%% patch size l1*l2
%l1=8;l2=8;
%

[n1,n2]=size(A);

if mode==1 %possible for other patching options

if nargin==2
   l1=8;l2=8;o1=4;o2=4; 
end

if nargin==4
   o1=round(l1/2);o2=round(l2/2);
end

tmp=mod(n1-l1,o1);
if tmp~=0
   A=[A;zeros(o1-tmp,n2)]; 
end
tmp=mod(n2-l2,o2);
if tmp~=0
   A=[A,zeros(size(A,1),o2-tmp)]; 
end



[N1,N2]=size(A);
 X=[];
for i1=1:o1:N1-l1+1
    for i2=1:o2:N2-l2+1
        tmp=reshape(A(i1:i1+l1-1,i2:i2+l2-1),l1*l2,1);
        X=[X,tmp];  
    end
end   
    
end


end

