function [Dinit]=yc_initD(l,c)
%yc_init: initialize ditionary
% BY Yangkang Chen
% Feb, 2021
%
% INPUT
%   din: input data
%   mode: patching mode
%   l: [l1,l2,l3]
%   l1: first patch size
%   l2: second patch size
%   l3: third patch size
%   s: [s1,s2,s3]
%   s1: first shifting size
%   s2: second shifting size
%   s3: third shifting size
%
%
% OUTPUT
% dout:
% Dinit: initial dictionary
%
l1=l(1);
l2=l(2);
l3=l(3);

c1=c(1);
c2=c(2);
c3=c(3);

%initialization
%[c1,c2,c3]: redundancy of the initial atom in 1st,2nd,3rd dimensions
%[l1,l2,l3]: patch sizes and the atom sizes in each dimension

dct1=zeros(l1,c1);
for k=0:1:c1-1
    V=cos([0:1:l1-1]'*k*pi/c1);
    if k>0
        V=V-mean(V);
    end
    dct1(:,k+1)=V/norm(V);
end

dct2=zeros(l2,c2);
for k=0:1:c2-1
    V=cos([0:1:l2-1]'*k*pi/c2);
    if k>0
        V=V-mean(V);
    end
    dct2(:,k+1)=V/norm(V);
end

if l3~=1
    dct3=zeros(l3,c3);
    for k=0:1:c3-1
        V=cos([0:1:l3-1]'*k*pi/c3);
        if k>0
            V=V-mean(V);
        end
        dct3(:,k+1)=V/norm(V);
    end
end

if l3==1
    Dinit=kron(dct1,dct2);%2D Dinit dictionary (l1*l2,c1*c2)
else
    Dinit=kron(kron(dct1,dct2),dct3);%3D Dinit dictionary (l1*l2*l3,c1*c2*c3)
end

return