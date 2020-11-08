function [dout]=yc_ddtf_recon(d0,mask,mode,l,s,perc,niter,a,param)
%yc_ddtf_recon: DDTF algorithm for 2D/3D reconstruction
% BY Yangkang Chen and Hang Wang
% Jan, 2020
%
% INPUT
%   mode: patching mode
%   l: [l1,l2,l3]
%   l1: first patch size
%   l2: second patch size
%   l3: third patch size
%   s: [s1,s2,s3]
%   s1: first shifting size
%   s2: second shifting size
%   s3: third shifting size
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
%
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, % Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.

[n1,n2,n3,n4,n5]=size(d0);

if isfield(param,'d0')
    d1=param.d0;
    fprintf('Using given model as the initilization\n');
else
    if isfield(param,'init') && (param.init==1)
        d1=yc_dlrecon_init(d0,mask,'nearest');
    fprintf('Using linear interpolation as the initilization\n');
    else
    d1=d0;
    fprintf('Using the input as the initilization\n');
    end
end


l1=l(1);
l2=l(2);
l3=l(3);

s1=s(1);
s2=s(2);
s3=s(3);

%initialization
%[c1,c2,c3]: redundancy of the initial atom in 1st,2nd,3rd dimensions
%[l1,l2,l3]: patch sizes and the atom sizes in each dimension
if ~isfield(param,'D')
    if isfield(param,'K')
        if n3==1
            c1=ceil(sqrt(param.K));
            c2=c1;
        else
            c1=ceil(param.K.^(1/3.0));
            c2=c1;
            c3=c1;
        end
    else
        c1=l1;
        c2=l2;
        c3=l3;
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
    
    if n3~=1
        dct3=zeros(l3,c3);
        for k=0:1:c3-1
            V=cos([0:1:l3-1]'*k*pi/c3);
            if k>0
                V=V-mean(V);
            end
            dct3(:,k+1)=V/norm(V);
        end
    end
    
    if n3==1
        DCT=kron(dct1,dct2);%2D DCT dictionary (64,256)
    else
        DCT=kron(kron(dct1,dct2),dct3);%3D DCT dictionary (l1*l2*l3,256)
    end
    param.D=DCT;
end


for iter=1:niter
    if n3==1
        X=yc_patch(d1,mode,l1,l2,s1,s2);
        [Dsgk,Gsgk]=yc_ddtf(X,param);
        Gsgkc=Gsgk;
        Gsgk=yc_pthresh(Gsgkc,'ph',perc);
        X2=Dsgk*Gsgk;
        d2=yc_patch_inv(X2,mode,n1,n2,l1,l2,s1,s2);
    else
        X=yc_patch3d(d1,mode,l1,l2,l3,s1,s2,s3);
        [Dsgk,Gsgk]=yc_ddtf(X,param);
        Gsgkc=Gsgk;
        Gsgk=yc_pthresh(Gsgkc,'ph',perc);
        X2=Dsgk*Gsgk;
        d2=yc_patch3d_inv(X2,mode,n1,n2,n3,l1,l2,l3,s1,s2,s3);
    end
    d1=a(iter)*d0.*mask+(1-a(iter))*d2.*mask+d2.*(1-mask);
    fprintf('Iter=%d/%d is finished\n',iter,niter);
end
dout=d1;

return