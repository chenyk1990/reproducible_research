function k = yc_kurtosis(x,flag)
% yc_kurtosis: Kurtosis
% By Yangkang Chen
% Jan, 2020
% INPUT
% x:     input
% flag:  bias (1) or unbias (0)   
% OUTPUT
% k: Kurtosis
%
% DEMO
% 1. a=randn(10,1);norm(yc_kurtosis(a)-kurtosis(a))
% 2. a=randn(10,1);norm(yc_kurtosis(a,0)-kurtosis(a,0))
% 
% assuming x is a vector
if nargin==1
   flag=1; 
end

[n,N]=size(x);
m=mean(x);

k=mean((x-ones(n,1)*m).^4)./(mean((x-ones(n,1)*m).^2).^2);

if flag==1
   k=k; 
else
   k=(n-1)*((n+1)*k-3*(n-1))/(n-2)/(n-3)+3; 
end

return

