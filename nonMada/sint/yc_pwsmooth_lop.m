function [dout] = yc_pwsmooth_lop(din,par,adj,add )
% yc_pwsmooth_lop: plane-wave smoothing operator
% BY Yangkang Chen, Aug, 10, 2021
% INPUT
% d: model/data
% par: parameter (dip,w1,ns,order,eps,nm,nd)
% adj: adj flag
% add: add flag
% OUTPUT
% m: data/model
%
% dottest:
% nt=501;nx=256;dip=randn(nt,nx);ns=10;order=1;eps=0.01;
% par=struct;par.dip=dip;par.w1=pwsmooth_set(dip,nt,nx,ns,order,eps);par.ns=ns;par.order=order;par.eps=eps;par.nm=nt*nx;par.nd=nt*nx;
% m1=randn(nt*nx,1);d1=[];
% [d1]=yc_pwsmooth_lop(m1,par,0,0);
% d2=randn(nt*nx,1);m2=[];
% m2=yc_pwsmooth_lop(d2,par,1,0);
% dot1=sum(sum(d1.*d2))
% dot2=sum(sum(m1.*m2))
%
% DEMO
% dn=levents(4);d0=fxymssa(dn,0,120,0.004,4,0);[nt,nx]=size(dn);dip=yc_dip2d_i(d0);ns=4;order=1;eps=0.01;
% par=struct;par.dip=dip;par.w1=pwsmooth_set(dip,nt,nx,ns,order,eps);par.ns=ns;par.order=order;par.eps=eps;par.nm=nt*nx;par.nd=nt*nx;
% [d1]=yc_pwsmooth_lop(dn,par,0,0);d1=reshape(d1,nt,nx);
% figure;yc_imagesc([dn,d1,dn-d1]);
%
% mada/test_pwsmooth.m

dip=par.dip;
w1=par.w1;
[n1,n2]=size(dip);
ns=par.ns;
order=par.order;
eps=par.eps;
ndn=n1*n2;
nds=n1*n2;
par.nm=ndn;
par.nd=nds;

if adj==1
    d=din;
    if isfield(par,'m') && add==1
        m=par.m;
    else
        m=zeros(par.nm,1);
    end
else
    m=din;
    if isfield(par,'d') && add==1
        d=par.d;
    else
        d=zeros(par.nd,1);
    end
end

% if adj
%     dn=m;
%     ds=d;
% else
%     dn=m;
%     ds=d;
% end

dn=m;
ds=d;

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
% fprintf("in pwsmooth_lop, size(w1)=(%d,%d)\n",size(w1,1),size(w1,2));
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

if adj==1
    dout=dn;
else
    dout=ds;
end

return



