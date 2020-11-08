function [D,G]=yc_sgk(X,param)
% yc_ksvd: K-SVD algorithm
% BY Yangkang Chen
% Jan, 2020
%
% INPUT
% X:     input training samples
% param: parameter struct
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
% DEMO: test/test_yc_sgk.m
%
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, % Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.

T=param.T;%T=1; %requred by SGK
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
        G=yc_ompN(D,X,T);
    else
        %error defined sparse coding
    end
   
    for ik=1:K %KSVD iteration, K times SVD
        inds=find(G(ik,:)~=0);
        if ~isempty(inds)
            D(:,ik)=sum(X(:,inds),2);
            D(:,ik)=D(:,ik)/norm(D(:,ik));
        end
    end
end

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

function [ G2 ] = yc_spray( G, n )
%multi-column sparse coding
% BY Yangkang Chen
% Jan, 2020
% G: row or column vector
% axis: 1 or 2
% n:  size

[n1,n2]=size(G);
if n1==1 %row vector, axis=1
    G2=zeros(n,n2);
    for i1=1:n
        G2(i1,:)=G;
    end
else %column vector, axis=2
    G2=zeros(n1,n);
    for i2=1:n
        G2(:,i2)=G;
    end
end
return
