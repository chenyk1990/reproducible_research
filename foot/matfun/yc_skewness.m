function s = yc_skewness(x,flag)
% yc_skewness: skewness
% By Yangkang Chen
% Jan, 2020
% INPUT
% x:     input
% flag:  bias (1) or unbias (0)   
% OUTPUT
% s: skewness
%
% DEMO
% 1. a=randn(10,1);norm(yc_skewness(a)-skewness(a))
% 2. a=randn(10,1);norm(yc_skewness(a,0)-skewness(a,0))
% 
% assuming x is a vector

if nargin==1
   flag=1; 
end
[n,N]=size(x);
m=mean(x);

s=mean((x-ones(n,1)*m).^3)./(sqrt(mean((x-ones(n,1)*m).^2)).^3);

if flag==1
   s=s; 
else
   s=sqrt(n*(n-1))/(n-2)*s; 
end

return

