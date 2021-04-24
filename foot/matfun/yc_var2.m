function k = yc_var2(D,N1,N2,axis,flag)
% yc_var2: variation calculation for a 2D Dictionary 
% By Wei Chen and Yangkang Chen
% Jan, 2021
% INPUT
% x:     input
% flag:  bias (1) or unbias (0)   
% OUTPUT
% k: Kurtosis
%
% DEMO
% 1. a=randn(10,1);norm(yc_kurtosis(a)-kurtosis(a))
% 2. a=randn(10,1);norm(yc_kurtosis(a,0)-kurtosis(a,0))
% 
% assuming x is a vector
if nargin==3
   axis=1;
   flag=1;
end

[n1,n2]=size(D);
if n1~=N1*N2
   error('dimension mismatch'); 
end

k=zeros(1,n2);
for i2=1:n2
 
    dtmp=reshape(D(:,i2),N1,N2);
    if axis==1
       k(i2)=sum(yc_var(dtmp).^2); 
    else
       k(i2)=sum(yc_var(dtmp').^2);  
    end
end


return

