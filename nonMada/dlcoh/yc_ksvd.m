function [D,G]=yc_ksvd(X,param)
%yc_ksvd: K-SVD algorithm
% BY Yangkang Chen
% Jan, 2020
% 
% INPUT
% X:     input training samples
%   param: parameter struct
%   param.mode=1;   %1: sparsity; 0: error
%   param.niter=10; %number of K-SVD iterations to perform; default: 10
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
% DEMO: test/test_yc_ksvd.m
% 
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, 
% Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.

T=param.T;
niter=param.niter;
mode=param.mode;
if isfield(param,'K')
K=param.K; 
else
K=size(D,2);%dictionary size: number of atoms
end
D=param.D(:,1:K);

for iter=1:niter
    
    if mode==1
        G=yc_ompN(D,X,T );
    else
        %error defined sparse coding
    end
    
    E0=X-D*G;%error before updating
    for ik=1:K %KSVD iteration, K times SVD
        %         inds=1:K;inds(ik)=[];
        E=E0+D(:,ik)*G(ik,:);
        inds=find(G(ik,:)~=0);
        R=E(:,inds);
        [u,s,v]=svds(R,1);
%         [U,S,V]=svd(R);
%         u=U(:,1);s=S(1,1);v=V(:,1);
        if ~isempty(u)
            D(:,ik)=u;
            G(ik,inds)=s*v.';
        end
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
