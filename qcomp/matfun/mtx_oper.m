function [ out ] = mtx_oper( in, Param, operator )
% General form of matrix based operator
% 
% By Yangkang Chen, Feb, 2015
% The University of Texas at Austin
% 
% operator = 1 : m->d (W )
% operator = -1: d->m (W')

W=Param.W;

if operator == 1
    out=W*in;
end

if operator == -1
    out=W'*in;
end
    

end

