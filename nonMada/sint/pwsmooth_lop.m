function [dn,ds] = pwsmooth_lop(adj,add,dip,w1,n1,n2,ns,order,eps,ndn,nds,dn,ds)
% pwsmooth_lop: Plane-wave smoothing (benchmarked with Madagascar, exact
% the same) 
% 
% This code is prepared with examples and dottest (exact)
% 
% Yangkang Chen, Zhejiang University, 2018
% 
% INPUT:
%
% adj:      adjoint flag
% add:      adding flag
% dip: slope (2D array)
% w1:       weight
% n1:       trace length
% n2:       number of traces
% ns:       spray radius
% order:    PWD order
% eps: regularization (default:0.01);
% 
% ndn: size of dn (n1*n2)
% nds: size of ds (n1*n2)
% 
% dn: model  (1D array) noisy data
% ds: data   (1D array)smoothed data
% 
% OUTPUT:
%
% dn: model
% ds: data
%
%  ALSO SEE
%  ~/chenyk/matlibcyk/test/test_structure_filtering.m
% 
% Example:
% d=fevents(2);[n1,n2]=size(d);figure;imagesc(d);
% dip=zeros(size(d));
% adj=0;add=0;[nt,nx]=size(d);ns=5;n=nt*nx;nu=nt*nx*(2*nr+1);order=1;eps=0.01;
% 
% [w1] = pwsmooth_set(dip,n1,n2,ns,order,eps);
% ndn=n1*n2;nds=n1*n2;
% [d,ds] = pwsmooth_lop(0,0,dip,w1,n1,n2,ns,order,eps,ndn,nds,d,ds);
% figure;imagesc(reshape(ds,n1,n2));
% 
% DOTTEST:
% DOT PRODUCT TEST TO PROVE THAT the operators are like a pair
% A and A'
% d=fevents(200);[n1,n2]=size(d);figure;imagesc(d);
% dip=zeros(size(d));
% adj=0;add=0;[nt,nx]=size(d);ns=5;n=nt*nx;nu=nt*nx*(2*nr+1);order=1;eps=0.01;
% [w1] = pwsmooth_set(dip,n1,n2,ns,order,eps);
% ndn=n1*n2;nds=n1*n2;ds=[];
% [d,ds] = pwsmooth_lop(0,0,dip,w1,n1,n2,ns,order,eps,ndn,nds,d,ds);
% figure;imagesc(reshape(ds,n1,n2));
% m1=randn(nt*nx,1);d1=[];
% [m11,d1]=pwsmooth_lop(0,0,dip,w1,n1,n2,ns,order,eps,n1*n2,n1*n2,m1,d1);%ACUTALLY M1=M11
% d2=randn(nt*nx,1);m2=[];
% [m2,d22]=pwsmooth_lop(1,0,dip,w1,n1,n2,ns,order,eps,n1*n2,n1*n2,m2,d2);%ACUTALLY M1=M11
% dot1=sum(sum(d1.*d2))
% dot2=sum(sum(m1.*m2))


if ndn~=nds
    error('Wrong size %d != %d',ndn,nds);
end

[ dn,ds ] = yc_adjnull( adj,add,ndn,nds,dn,ds );


ns2=2*ns+1;%spray diameter
n12=n1*n2;

u=zeros(n1,ns2,n2);
utmp=zeros(n12*ns2,1);
w=zeros(ns2,1);
% w1=zeros(n1,n2);

for is=0:ns2-1
    w(is+1)=ns+1-abs(is-ns);
end


% for Normalization
% t=zeros(n12,1);

if adj

    for i2=0:n2-1
       for i1=0:n1-1
            ws=w1(i1+1,i2+1);
            for is=0:ns2-1
                u(i1+1,is+1,i2+1)=ds(i2*n1+i1+1)*w(is+1)*ws;

            end
       end
    end  

    
    utmp=u(:);
    

    [dn,utmp]=pwspray_lop(1,1,n12,n12*ns2,dn,utmp,dip,ns,n1,n2,order,eps);
    

else

    [dn,utmp]=pwspray_lop(0,0,n12,n12*ns2,dn,utmp,dip,ns,n1,n2,order,eps);
 
    u=reshape(utmp,n1,ns2,n2);
    
    for i2=0:n2-1
        for i1=0:n1-1
            ws=w1(i1+1,i2+1);
            for is=0:ns2-1
                ds(i2*n1+i1+1)=ds(i2*n1+i1+1)+u(i1+1,is+1,i2+1)*w(is+1)*ws;
            end
        end
    end
end


return


%






