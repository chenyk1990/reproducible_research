function [D,G]=yc_sgk(X,param)
%yc_sgk: SGK algorithm
% BY Yangkang Chen
% Jan, 2020
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
% DEMO: test/test_yc_sgk.m
% 
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, 
% Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.

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
        G=ompN(D,X,1);
        % exact form
        % 
        
    else
        %error defined sparse coding
    end
    
    for ik=1:K %KSVD iteration, K times SVD
        inds=find(G(ik,:)~=0);
        if ~isempty(inds)
            D(:,ik)=sum(X(:,inds),2);%better using a weighted summation ? NO, equivalent 
            D(:,ik)=D(:,ik)/norm(D(:,ik)); 
        end
    end
end

%% extra step
G=ompN(D,X,T);

return

function [ G ] = ompN( D, X, K )
%multi-column sparse coding
% BY Yangkang Chen
% Jan, 2020
[n1,n2]=size(X);
[n1,n3]=size(D);
G=zeros(n3,n2);
if K==1
    for i2=1:n2
        G(:,i2)=omp_e(D,X(:,i2));
    end
else
    for i2=1:n2
        G(:,i2)=yc_omp0(D,X(:,i2),K);
    end
end
return

function [ g ] = omp_e( D, x )
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
g=zeros(n2,1);

k=1;
max=0;                      %initialize max=0
for i2=1:n2                 %loop over all atoms (greedy algorithm)
    dtr=abs(sum(D(:,i2).*x));
    if max<dtr
        max=dtr;
        k=i2;
    end
end

g(k)=sum(D(:,k).*x)/sum(D(:,k).*D(:,k)); %option 1: more accurate sparse coding (according to definition of sparse coding)
% g(k)=1.0; %option 2: according to the requirement of equation (8) in Chen (2017), GJI.
% Note: option1 and option2 result in exactly the same result. 
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
