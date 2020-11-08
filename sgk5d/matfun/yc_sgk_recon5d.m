function [dout]=yc_sgk_recon5d(d0,mask,mode,l,s,perc,niter,a,param)
%yc_sgk_recon5d: SGK algorithm for 4D/5D reconstruction
% BY Yangkang Chen
% April, 2020 (3D-4D/5D)
%
% INPUT
%   mode: patching mode
%   l: [l1,l2,l3,l4,l5]
%   l1: first patch size
%   l2: second patch size
%   l3: third patch size
%   l4: fourth patch size
%   l5: fifth patch size
%   s: [s1,s2,s3,s4,s5]
%   s1: first shifting size
%   s2: second shifting size
%   s3: third shifting size
%   s4: fourth shifting size
%   s5: fifth shifting size
%   perc: percentage
% param: parameter struct for DL
%   param.mode=1;   %1: sparsity; 0: error
%   param.niter=10; %number of SGK iterations to perform; default: 10
%   param.D=DCT;    %initial D
%   param.T=3;      %sparsity level
%
% OUTPUT
% dout:
%
% for X=DG
% size of X: MxN
% size of D: MxK
% size of G: KxN
%
% DEMO:
% test/test_yc_sgk_recon.m
% test/test_yc_sgk_recon3d.m
% test/test_yc_sgk_recon4d.m
% test/test_yc_sgk_recon5d.m
% 
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, % Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.
[n1,n2,n3,n4,n5]=size(d0);
d1=d0;

l1=l(1);
l2=l(2);
l3=l(3);
l4=l(4);
l5=l(5);

s1=s(1);
s2=s(2);
s3=s(3);
s4=s(4);
s5=s(5);

%initialization
%[c1,c2,c3,c4,c5]: redundancy of the initial atom in 1st,2nd,3rd,4th,5th dimensions
%[l1,l2,l3,l4,l5]: patch sizes and the atom sizes in each dimension
if ~isfield(param,'D')
    if isfield(param,'K')
        if n5==1
            c1=ceil(param.K.^(1/4.0));
            c2=c1;
            c3=c1;
            c4=c1;
        else
            c1=ceil(param.K.^(1/5.0));
            c2=c1;
            c3=c1;
            c4=c1;
            c5=c1;
        end
    else
        c1=l1;
        c2=l2;
        c3=l3;
        c4=l4;
        c5=l5; 
    end
    
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

    dct3=zeros(l3,c3);
    for k=0:1:c3-1
        V=cos([0:1:l3-1]'*k*pi/c3);
        if k>0
            V=V-mean(V);
        end
        dct3(:,k+1)=V/norm(V);
    end
    
    dct4=zeros(l4,c4);
    for k=0:1:c4-1
        V=cos([0:1:l4-1]'*k*pi/c4);
        if k>0
            V=V-mean(V);
        end
        dct4(:,k+1)=V/norm(V);
    end
    
    if n5~=1
        dct5=zeros(l5,c5);
        for k=0:1:c5-1
            V=cos([0:1:l5-1]'*k*pi/c5);
            if k>0
                V=V-mean(V);
            end
            dct5(:,k+1)=V/norm(V);
        end
    end
    
    if n5==1
        DCT=kron(kron(kron(dct1,dct2),dct3),dct4);
    else
        DCT=kron(kron(kron(kron(dct1,dct2),dct3),dct4),dct5);
    end
    param.D=DCT;
end


for iter=1:niter
    X=yc_patch5d(d1,mode,l1,l2,l3,l4,l5,s1,s2,s3,s4,s5);
    [Dsgk,Gsgk]=yc_sgk(X,param);
    Gsgkc=Gsgk;
    Gsgk=yc_pthresh(Gsgkc,'ph',perc);
    X2=Dsgk*Gsgk;
    d2=yc_patch5d_inv(X2,mode,n1,n2,n3,n4,n5,l1,l2,l3,l4,l5,s1,s2,s3,s4,s5);
    d1=a(iter)*d0.*mask+(1-a(iter))*d2.*mask+d2.*(1-mask);
    fprintf('Iter=%d/%d is finished\n',iter,niter);
end
dout=d1;

return
