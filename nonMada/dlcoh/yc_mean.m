function m = yc_mean(x)
% yc_mean: mean
% By Yangkang Chen
% Jan, 2020
% INPUT
% x:     input
% flag:  bias (1) or unbias (0)   
% OUTPUT
% k: mean
%
% DEMO
% 1. a=randn(10,1);norm(yc_mean(a)-mean(a))
% 
% assuming x is a vector

m=mean(x);

return

