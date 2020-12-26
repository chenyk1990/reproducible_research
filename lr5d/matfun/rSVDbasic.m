function [dout]=rSVDbasic(A, k, P, F)
%function [U, S, V]=rSVDbasic(A, k, P)
%Basic randomized truncated singular value decomposition of A.
%Syntax:
% s = rSVDbasic(A, k)
% s = rSVDbasic(A, k, P)
% [U, S, V]= rSVDbasic(A, k)
% [U, S, V]= rSVDbasic(A, k, P)
% -P is an optional parameter to balance time and accuacy (default value 0).
%  With large P, the accuracy increases with runtime overhead.
%Algorithm: the basic randomized scheme for computing truncated SVD.
%  It outputs the largest k singular values, or the corresponding factors.

if nargin<3,
    P=0;
end

s=10;               % over-sampling
[m,n]= size(A);
B= randn(n, k+s);
U= A*B;
[U, ~]= qr(U, 0);
for j=1:P,
    [B, ~]= qr(A'*U, 0);   % May reduce an orthogonalization
    [U, ~]= qr(A*B, 0);
end
B= A'*U;

% if nargout==1,
%     U= svd(B','econ');
%     U= U(1:k);
% else
    [U1, S, V]= svd(B', 'econ');
    U= U*U1(:,1:k);
%     S= S(1:k,1:k);
    
    for j=1:k
        S(j,j)=S(j,j)*(1-S(k+1,k+1)^F/S(j,j)^F);
    end 
    
    S= S(1:k,1:k);
%     S= diag(S);
%     S= S(1:k);
    V= V(:,1:k);
% end

dout = U*S*V';

end 