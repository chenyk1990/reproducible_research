function [w,diag,offd,trace] = predict_step(e,nw,adj,forw,n1,pp,trace)
%predict_step: prediction ste
%
% INPUT:
%
% e: regularization parameter (default, 0.01*0.01);
% nw: accuracy order
% two: if two predictions (neglected)
% adj: adjoint flag
% forw: forward or backward
% n1: trace length
% pp: slope
% trace: input trace
% 
% OUTPUT:
% w: PWD object
% diag,offd: diagonal/offdiagonals of the banded matrix
% trace: output trace

nb=2*nw;
eps=e;
eps2=e;

diag=zeros(n1,1);
offd=zeros(n1,nb);


[diag,offd] = regularization(diag,offd,nw,eps,eps2);


% b = banded_solve(n,band,diag,offd,b)

%only define diagonal,offdiagonal

[w,diag,offd] = pwd_define(forw,diag,offd,n1,nw,pp);
% [diag,offd]
   
if adj
   trace = banded_solve(n1,nb,diag,offd,trace);
end

t0=trace(1);
t1=trace(2);
t2=trace(n1-1);
t3=trace(n1);


[trace] = pwd_set(adj,w,diag,offd,pp,trace);%?


trace(1)=trace(1)+eps2*t0;
trace(2)=trace(2)+eps2*t1;

trace(n1-1)=trace(n1-1)+eps2*t2;
trace(n1)=trace(n1)+eps2*t3;

if ~adj
   trace = banded_solve(n1,nb,diag,offd,trace);   
end



return


