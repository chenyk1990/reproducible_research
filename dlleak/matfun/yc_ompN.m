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