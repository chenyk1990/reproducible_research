function [D,G]=yc_ddtf(X,param)
%yc_ddtf: DDTF algorithm
% BY Yangkang Chen
% Jan, 2020
% Modified April, 2020
%
% INPUT
% X:     input training samples
% param: parameter struct
%   param.mode=1;   %1: sparsity; 0: error
%   param.niter=10; %number of SGK iterations to perform; default: 10
%   param.D=DCT;    %initial D
%   param.T=3;      %sparsity level
%
% OUTPUT
% D:    learned dictionary
% G:    sparse coefficients
%
% for X=DG
% size of X: MxN
% size of D: MxK
% size of G: KxN
%
% SEE ALSO: yc_ddtf0.m (initial version)
% DEMO: test/test_yc_ddtf_denoise.m
%
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, 
% Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.

M=size(X,1);
T=param.T;
niter=param.niter;
mode=param.mode;
if isfield(param,'K')
K=param.K; 
else
K=size(D,2);%dictionary size: number of atoms
end
D=param.D(:,1:K);

if isfield(param,'lambda')
    lambda=param.lambda;
else
    lambda=0.5;
end

for iter=1:niter
    
    if mode==1
        G=yc_ompN(D,X,T );
    else
        %error defined sparse coding
        if mode==0
        G = D'*X;
        G = yc_wthresh(G, 'h', lambda);    
        end
    end


    XG = X*G';           %size of X is M*N; size of G is K*N
    [U, ~, V] = svd(XG); %size of U is M*M; size of V is K*K
    if M<=K
    D = U*V(:,1:M)';
    else
    D = U(:,1:K)*V';    
    end
end

%% extra step
G=yc_ompN(D,X,T);

return

function [ G ] = yc_ompN( D, X, K )
%multi-column sparse coding
% BY Yangkang Chen
% Jan, 2020
[n1,n2]=size(X);
[n1,n3]=size(D);
G=zeros(n3,n2);
for i2=1:n2
    G(:,i2)=yc_omp0(D,X(:,i2),K);
end

return

function [ g ] = yc_omp0( D, x, K )
% yc_omp0: Most basic orthogonal matching pursuit for sparse coding
% this simple tutorial code use the sparseity-constrained sparse coding
% model
% 
% two version of sparse coding problems:
% 1) Sparseity-constrained sparse coding problem
%   gamma = min |x-Dg|_2^2  s.t. |g|_0 <= K 
% 2) Error-constrained sparse coding problem
%   gamma = min |g|_0       s.t. |x-Dg|_2^2 <=eps
% 
% Author: Yangkang Chen
% Oct 25, 2016
[n1,n2]=size(D);
I=[];
r=x;
g=zeros(n2,1);
for ik=1:K
    k=[];
    max=0;                      %initialize max=0
    for i2=1:n2                 %loop over all atoms (greedy algorithm)       
        if sum(find(I==i2))==0  %search among the other atoms
            dtr=abs(sum(D(:,i2).*r));
            if max<dtr
                max=dtr;
                k=i2;
            end        
        end
    end
    I=[I,k];
    g(I)=inv(D(:,I)'*D(:,I))*D(:,I)'*x; %g_I = D_I^{+}x, D_I^TD_I is guaranteed to be non-singular 
    r=x-D(:,I)*g(I);
end

return
