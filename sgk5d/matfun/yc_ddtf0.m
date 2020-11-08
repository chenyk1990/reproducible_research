function [D,G]=yc_ddtf0(X,param)
%yc_ddtf0: DDTF algorithm
% BY Yangkang Chen and Hang Wang
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
% DEMO: test/test_yc_ddtf.m
%
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, % Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.

T=param.T;
niter=param.niter;
mode=param.mode;
if isfield(param,'K')
    K=param.K;
else
    K=size(D,2);%dictionary size: number of atoms
end


if isfield(param,'lambda')
    lambda=param.lambda;
else
    lambda=0.5;
end


if isfield(param,'D') && ~isempty(param.D)
    D=param.D(:,1:K);
else
    
    %% initialize
    if isfield(param,'htype')
        htype=param.htype;
    else
        htype='spline';
    end
    
    if isfield(param,'level')
        level=param.level;
    else
        level=2;
    end
    switch htype
        
        case 'spline' % Piece-wise Linear Spline
            [h, hsize] = TFDict(level);
            
        case 'haar' % Haar transform matrix
            hsize = 2^level;
            h = haarmtx(hsize);
            
        case 'dct' % DCT transform matrix
            hsize = 2^level;
            h = dctmtx(hsize);
    end
    
    D = kron(h', h');
end

for iter=1:niter
    
    
    %     if mode==1
    %         G=yc_ompN(D,X,T );
    %     else
    %         %error defined sparse coding
    %     end
    
    G = D'*X;
    G = yc_wthresh(G, 'h', lambda);
    
    XG = X*G';
    [U, ~, V] = svd(XG);
    
    D = U*V';
    
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


%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H ] = haarmtx(n)
% k = 2^p + q
% p, specifies the magnitude and width (or scale) of the shape;
% q, specifies the position (or shift) of the shape.

H = zeros(n,n);
H(1,1:n) = ones(1,n)/sqrt(n);

for k = 1:n-1
    
    p = fix(log(k)/log(2));
    
    q = k-(2^p);
    
    k1 = 2^p;
    t1 = n/k1;
    
    k2 = 2^(p+1);
    t2 = n/k2;
    
    for i = 1:t2
        H(k+1,i+q*t1) = (2^(p/2))/sqrt(n);
        H(k+1,i+q*t1+t2) = -(2^(p/2))/sqrt(n);
    end
    
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H, hsize] = TFDict(level)
% Generate Tight Frame Dictionary
% input:
%  level: level of Tight Frame
% output:
%  dict: Tight Frame dictionary

if nargin < 1, level = 1; end

phi = 1/4*[1,2 ,1]; % refinement mask
psi1 = sqrt(2)/4*[1,0,-1]; % first wavelet mask
psi2 = 1/4*[-1,2,-1]; % second wavelet mask

a0 = phi';
a1 = psi1';
a2 = psi2';
if level == 1
    hsize = 3;
    H = [phi;psi1;psi2];
elseif level == 2
    h0 = conv([a0(1),0,a0(2),0,a0(3)],a0);
    h1 = conv([a1(1),0,a1(2),0,a1(3)],a0);
    h2 = conv([a2(1),0,a2(2),0,a2(3)],a0);
    h3 = [a1',0,0,0,0];
    h4 = [0,0,0,0,a1'];
    h5 = [a2',0,0,0,0];
    h6 = [0,0,0,0,a2'];
    
    hsize = 7;
    H = [h0;h1;h2;h3;h4;h5;h6];
elseif level == 3
    h0 = conv([a0(1),0,a0(2),0,a0(3)],a0);
    h1 = conv([a1(1),0,a1(2),0,a1(3)],a0);
    h2 = conv([a2(1),0,a2(2),0,a2(3)],a0);
    h3 = [a1',0,0,0,0]*sqrt(2)/2;
    h4 = [0,0,0,0,a1']*sqrt(2)/2;
    h5 = [a2',0,0,0,0]*sqrt(2)/2;
    h6 = [0,0,0,0,a2']*sqrt(2)/2;
    
    I0 = conv([h0(1),0,h0(2),0,h0(3),0,h0(4),0,h0(5),0,h0(6),0,h0(7)],a0);
    I1 = conv([h1(1),0,h1(2),0,h1(3),0,h1(4),0,h1(5),0,h1(6),0,h1(7)],a0);
    I2 = conv([h2(1),0,h2(2),0,h2(3),0,h2(4),0,h2(5),0,h2(6),0,h2(7)],a0);
    I3 = [h1,0,0,0,0,0,0,0,0]*sqrt(2)/2;
    I4 = [0,0,0,0,0,0,0,0,h1]*sqrt(2)/2;
    I5 = [h2,0,0,0,0,0,0,0,0]*sqrt(2)/2;
    I6 = [0,0,0,0,0,0,0,0,h2]*sqrt(2)/2;
    I7 = [h3,0,0,0,0,0,0,0,0]*sqrt(2)/2;
    I8 = [0,0,0,0,0,0,0,0,h3]*sqrt(2)/2;
    I9 = [h4,0,0,0,0,0,0,0,0]*sqrt(2)/2;
    I10 = [0,0,0,0,0,0,0,0,h4]*sqrt(2)/2;
    I11 = [h5,0,0,0,0,0,0,0,0]*sqrt(2)/2;
    I12 = [0,0,0,0,0,0,0,0,h5]*sqrt(2)/2;
    I13 = [h6,0,0,0,0,0,0,0,0]*sqrt(2)/2;
    I14 = [0,0,0,0,0,0,0,0,h6]*sqrt(2)/2;
    
    hsize = 15;
    H = [I0;I1;I2;I3;I4;I5;I6;I7;I8;I9;I10;I11;I12;I13;I14];
end

% dict = kron(H,H);
return


