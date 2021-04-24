function v = yc_var(x)
% yc_var: var
% By Yangkang Chen
% Jan, 2020
% INPUT
% x:     input
% flag:  bias (1) or unbias (0)   
% OUTPUT
% k: var
%
% DEMO
% 1. a=randn(10,1);norm(yc_var(a)-var(a))
% 
% assuming x is a vector
n=size(x,1);

m=mean(x);

v=sum((x-m).^2)/(n-1); %% very important (n-1) 

return