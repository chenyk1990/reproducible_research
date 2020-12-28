function m = yc_median(x)
% yc_median: median
% By Yangkang Chen
% Jan, 2020
% INPUT
% x:     input
% flag:  bias (1) or unbias (0)   
% OUTPUT
% m: median
%
% DEMO
% 1. a=randn(10,1);norm(yc_median(a)-median(a))
% 
% assuming x is a vector

m=median(x);

return

