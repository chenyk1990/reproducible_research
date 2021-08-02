function [diag,offd] = regularization(diag,offd,nw,eps,eps2)
% fill diag and offd using regularization */
%
% INPUT:
%
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% eps: regularization parameter (default: e*e, e=0.01);
% eps2: second regularization parameter (default, same as eps)
% nw: accuracy order (nb=2*nw)
%
% OUTPUT:
% 
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
%
nb=2*nw;
n1=length(diag);

for i1=0:n1-1
    diag(i1+1)=6*eps;
    offd(i1+1,1)=-4*eps;
    offd(i1+1,2)=eps;
    
    for ib=2:nb-1
        offd(i1+1,ib+1)=0.0;
    end
end

    diag(1)=eps2+eps;
    diag(2)=eps2+5*eps;
    
    diag(n1)=eps2+eps;
    diag(n1-1)=eps2+5*eps;  
    
    offd(1,1)=-2*eps;
    offd(n1-1,1)=-2*eps;

end


