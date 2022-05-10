function [u1,u] = pwspray_lop(adj,add,n,nu,u1,u,dip,nr,nt,nx,order,eps)
%pwspray_lop: Plane-wave spray operator(benchmarked with Madagascar, exact
% the same)
%
% This code is prepared with examples and dottest (exact)
%
%
% Yangkang Chen, Zhejiang University, 2018
%
% INPUT:
%
% adj: adjoint flag
% add: adding flag
% n: size of u1 (n = n1*n2, size of the input data matrix)
% nu: size of u (nu = n*(2*r+1), r is spray radius)
% u1: model  (1D array)
% u: data    (1D array)
% dip: slope (2D array)
% nr: spray radius
% nx: number of traces
% order: PWD order
% eps: regularization (default:0.01);
%
% OUTPUT:
%
% u1: model
% u: data
%
%
% DEMO1:
% d=fevents(200);figure;imagesc(d);
% dip=zeros(size(d));
% adj=0;add=0;[nt,nx]=size(d);nr=5;n=nt*nx;nu=nt*nx*(2*nr+1);order=1;eps=0.01;
% u1=d(:);u=[];
% [u1,u]=pwspray_lop(adj,add,n,nu,u1,u,dip,nr,nt,nx,order,eps);
%
% figure;imagesc(reshape(u1,nt,nx));
% u3d=reshape(u,nt,2*nr+1,nx);
% figure;imagesc(u3d(:,:,1));
% figure;imagesc(u3d(:,:,50));
%
% DEMO2:
% d=levents(200);figure;imagesc(d);
% [dip] = dip2d_shan(d,3,10,1,0.1,[10,10,1]);figure;imagesc(dip);colorbar;
% adj=0;add=0;[nt,nx]=size(d);nr=5;n=nt*nx;nu=nt*nx*(2*nr+1);order=1;eps=0.01;
% u1=d(:);u=[];
% [u1,u]=pwspray_lop(adj,add,n,nu,u1,u,dip,nr,nt,nx,order,eps);
%
% figure;imagesc(reshape(u1,nt,nx));
% u3d=reshape(u,nt,2*nr+1,nx);
%
% for i=1:50
%     figure(10);pause(0.1);imagesc(u3d(:,:,i));
% end
%
% DOTTEST:
% DOT PRODUCT TEST TO PROVE THAT the operators are like a pair
% A and A'
%
% d=fevents(200);figure;imagesc(d);
% dip=zeros(size(d));
% adj=0;add=0;[nt,nx]=size(d);nr=5;n=nt*nx;nu=nt*nx*(2*nr+1);order=1;eps=0.01;
% u1=d(:);u=[];
%
% m1=randn(nt*nx,1);d1=[];
% [m11,d1]=pwspray_lop(0,add,n,nu,m1,d1,dip,nr,nt,nx,order,eps);%ACUTALLY M1=M11
%
% d2=randn(nt*nx*(2*nr+1),1);m2=[];
% [m2,d22]=pwspray_lop(1,add,n,nu,m2,d2,dip,nr,nt,nx,order,eps);%ACUTALLY M1=M11
%
% dot1=sum(sum(d1.*d2))
% dot2=sum(sum(m1.*m2))
%

% %


n1=nt;
n2=nx;
ns=nr;     %spray radius
ns2=2*ns+1;%spray diameter
p=dip;
trace=zeros(n1,1);
e=eps*eps;
nw=order;

if n~=n1*n2
    error('Wrong size %d != %d*%d',n,n1,n2);
end

if nu~=n*ns2
    error('Wrong size %d != %d*%d',nu,n,ns2);
end

[ u1,u ] = yc_adjnull( adj,add,n,nu,u1,u );

for i=0:n2-1
    
    if adj
        for i1=0:n1-1
            trace(i1+1)=0.0;
        end
        
        %predict forward
        for is=ns-1:-1:0
            ip=i+is+1;
            if ip>=n2
                continue;
            end
            j=ip*ns2+ns+is+1;
            for i1=0:n1-1
                trace(i1+1)=trace(i1+1)+u(j*n1+i1+1);
            end
            [w,diag,offd,trace] = predict_step(e,nw,1,1,n1,p(:,ip),trace);
        end
        
        for i1=0:n1-1
            u1(i*n1+i1+1)= u1(i*n1+i1+1)+trace(i1+1);
            trace(i1+1)=0.0;
        end
        % predict backward
        for is=ns-1:-1:0
            ip=i-is-1;
            if ip<0
                continue;
            end
            j=ip*ns2+ns-is-1;
            for i1=0:n1-1
                trace(i1+1)=trace(i1+1)+u(j*n1+i1+1);
            end
            [w,diag,offd,trace] = predict_step(e,nw,1,0,n1,p(:,ip+1),trace);
        end
        
        for i1=0:n1-1
            u1(i*n1+i1+1)=u1(i*n1+i1+1)+trace(i1+1);
            trace(i1+1)=u((i*ns2+ns)*n1+i1+1);
            u1(i*n1+i1+1)=u1(i*n1+i1+1)+trace(i1+1);
        end
        
    else
        for i1=0:n1-1
            trace(i1+1)=u1(i*n1+i1+1);
            u((i*ns2+ns)*n1+i1+1) = u((i*ns2+ns)*n1+i1+1)+trace(i1+1);
        end
        
        %predict forward
        for is=0:ns-1
            ip=i-is-1;
            if ip<0
                break;
            end
            j=ip*ns2+ns-is-1;
            [w,diag,offd,trace] = predict_step(e,nw,0,0,n1,p(:,ip+1),trace);
            for i1=0:n1-1
                u(j*n1+i1+1)=u(j*n1+i1+1)+trace(i1+1);
            end
        end
        
        for i1=0:n1-1
            trace(i1+1)=u1(i*n1+i1+1);
        end
        
        % predict backward
        for is=0:ns-1
            ip=i+is+1;
            if ip>=n2
                break;
            end
            j=ip*ns2+ns+is+1;
            [w,diag,offd,trace] = predict_step(e,nw,0,1,n1,p(:,ip),trace);
            for i1=0:n1-1
                u(j*n1+i1+1)= u(j*n1+i1+1)+trace(i1+1);
            end
        end
    end
    
end

return



function [w,diag,offd,trace] = predict_step(e,nw,adj,forw,n1,pp,trace)
%predict_step: prediction ste
%
% INPUT:
%
% e: regularization parameter (default, 0.01*0.01);
% nw: accuracy order
% two: if two predictions (neglected)
% adj: adjoint flag
% forw: forward or backward
% n1: trace length
% pp: slope
% trace: input trace
%
% OUTPUT:
% w: PWD object
% diag,offd: diagonal/offdiagonals of the banded matrix
% trace: output trace

nb=2*nw;
eps=e;
eps2=e;

diag=zeros(n1,1);
offd=zeros(n1,nb);


[diag,offd] = regularization(diag,offd,nw,eps,eps2);


% b = banded_solve(n,band,diag,offd,b)

%only define diagonal,offdiagonal

[w,diag,offd] = pwd_define(forw,diag,offd,n1,nw,pp);
% [diag,offd]

if adj
    trace = banded_solve(n1,nb,diag,offd,trace);
end

t0=trace(1);
t1=trace(2);
t2=trace(n1-1);
t3=trace(n1);


[trace] = pwd_set(adj,w,diag,offd,pp,trace);%?


trace(1)=trace(1)+eps2*t0;
trace(2)=trace(2)+eps2*t1;

trace(n1-1)=trace(n1-1)+eps2*t2;
trace(n1)=trace(n1)+eps2*t3;

if ~adj
    trace = banded_solve(n1,nb,diag,offd,trace);
end



return


function [diag,offd] = regularization(diag,offd,nw,eps,eps2)
% fill diag and offd using regularization 
%
% INPUT:
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% eps: regularization parameter (default: e*e, e=0.01);
% eps2: second regularization parameter (default, same as eps)
% nw: accuracy order (nb=2*nw)
%
% OUTPUT:
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

nb=2*nw;
n1=length(diag);

for i1=0:n1-1
    diag(i1+1)=6*eps;
    offd(i1+1,1)=-4*eps;
    offd(i1+1,2)=eps;
    
    for ib=2:nb-1
        offd(i1+1,ib+1)=0.0;
    end
end

    diag(1)=eps2+eps;
    diag(2)=eps2+5*eps;
    
    diag(n1)=eps2+eps;
    diag(n1-1)=eps2+5*eps;  
    
    offd(1,1)=-2*eps;
    offd(n1-1,1)=-2*eps;
return

function [w,diag,offd] = pwd_define(forw,diag,offd,n1,nw,pp)
%pwd_define: matrix multiplication
%
% INPUT:
%
% forw:forward or backward
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% n1: trace length
% nw: PWD filter(accuracy) order (default nw=1)
% pp: slope                  (1D array)
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% ifset: only define diagonal,offdiagonal or apply pwd_set


%define PWD object(struct)
w=struct;
w.n=n1;
w.na=2*nw+1;
w.a=zeros(n1,w.na);%2D array
w.b=zeros(w.na,1);%1D array

%apfilt_init(nw)
% nw=(w.na-1)/2;
n=w.n;
nb=2*nw;

%     for (i=0; i < n; i++) {
% 	passfilter (pp[i], w->b);


%slv=sf_banded_init(n1,nb);
% diag=zeros(n1,1);
% offd=zeros(n1,nb);

for i=0:n-1
    w.b=passfilter(pp(i+1),nw);
    for j=0:w.na-1
        if(forw)
            w.a(i+1,j+1)=w.b(w.na-j);
        else
            w.a(i+1,j+1)=w.b(j+1);
        end
        
    end
end

for i=0:n-1
    for j=0:w.na-1
        k=i+j-nw;
        if k>=nw && k<n-nw
            aj=w.a(k+1,j+1);
            diag(i+1)=diag(i+1)+aj*aj;
        end
    end
    for m=0:2*nw-1
        for j=m+1:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                aj=w.a(k+1,j+1);
                am=w.a(k+1,j-m);
                offd(i+1,m+1)=offd(i+1,m+1)+am*aj;
            end
        end
    end
    
end

return


function [out] = pwd_set(adj,w,diag,offd,pp,inp)
%pwd_set: matrix multiplication
%
% INPUT:
%
% adj:adjoint flag
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% n1: trace length
% nw: PWD filter(accuracy) order (default nw=1)
% pp: slope                  (1D array)
% diag: defined diagonal .   (1D array) (useless)
% offd: defined off-diagonal (2D array) (useless)
% inp: model
% out: data


%define PWD object(struct)
% w=struct;
% w.n=n1;
% w.na=2*nw+1;
% w.a=zeros(n1,w.na);%2D array
% w.b=zeros(w.na,1);%1D array

%apfilt_init(nw)
% nw=(w.na-1)/2;
n=w.n;
nw=(w.na-1)/2;

%% pwd_set
tmp=zeros(n,1);

if adj
    for i=0:n-1
        tmp(i+1)=0.0;
    end
    for i=0:n-1
        for j=0:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                tmp(k+1)=tmp(k+1)+w.a(k+1,j+1)*inp(i+1);
            end
            
        end
    end
    
    for i=0:n-1
        out(i+1)=0.0;
    end
    
    for i=nw:n-nw-1
        for j=0:w.na-1
            k=i+j-nw;
            out(k+1)=out(k+1)+w.a(i+1,j+1)*tmp(i+1);
        end
    end
else
    for i=0:n-1
        tmp(i+1)=0.0;
    end
    for i=nw:n-nw-1
        for j=0:w.na-1
            k=i+j-nw;
            tmp(i+1)=tmp(i+1)+w.a(i+1,j+1)*inp(k+1);
        end
    end
    for i=0:n-1
        out(i+1)=0.0;
        for j=0:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                out(i+1)=out(i+1)+w.a(k+1,j+1)*tmp(k+1);
            end
        end
    end
    
end


return


function [a,b] = passfilter(p,nw)
% passfilter: find filter coefficients£¨verfied)
% All-pass plane-wave destruction filter coefficients
%
% INPUT:
% p: slope
% a: output filter (n+1) (1D array)
%
n=nw*2;
b=zeros(n+1,1);
for k=0:n
    bk=1;
    for j=0:n-1
        if (j<n-k)
            bk=bk*(k+j+1.0)/(2*(2*j+1)*(j+1));
        else
            bk=bk*1.0/(2*(2*j+1));
        end
        
    end
    b(k+1)=bk;
end

for k=0:n
    ak=b(k+1);
    for j=0:n-1
        if j<n-k
            ak=ak*(n-j-p);
        else
            ak=ak*(p+j+1);
        end
    end
    a(k+1)=ak;
end
return


function b = banded_solve(n,band,diag,offd,b)
%banded_solve: Banded matrix solver
%
% Yangkang Chen, Zhejiang University, Oct, 2018
%
% INPUT:
%
% n:    matrix size
% band: band size
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
%
% OUTPUT:
% b: solution
%
% Test:
%  A=[1,0.5,0.1,0;0.5,2,0.3,0.05;0.1,0.3,3.0,0.2;0,0.05,0.2,4.0];
%  b=[1,2,1.5,1.0]';
%  x1=A\b;
%  diag=[1,2,3,4];offd=[[0.5;0.3;0.2;0],[0.1;0.05;0;0]];
%  x2=banded_solve(4,2,diag,offd,b);
%  norm(x1-x2)
%
% Reference:
% Fomel, Madagascar, api/c/banded.c

%define Band object(struct)
slv=struct;
slv.n=n;
slv.band=band;
slv.d=zeros(n,1);       %1D array
slv.o=zeros(n-1,band);  %2D array


% define the banded matrix
for k=0:slv.n-1
    t=diag(k+1);
    m1=min(k,slv.band);
    for m=0:m1-1
        t=t-slv.o(k-m,m+1)*slv.o(k-m,m+1)*slv.d(k-m);
    end
    slv.d(k+1)=t;
    n1=min(slv.n-k-1,slv.band);
    for n=0:n1-1
        t=offd(k+1,n+1);
        m1=min(k,slv.band-n-1);
        for m=0:m1-1
            t=t-slv.o(k-m,m+1)*slv.o(k-m,n+m+2)*slv.d(k-m);
        end
        slv.o(k+1,n+1)=t/slv.d(k+1);
    end
    
    
end

% the solver
for k=1:slv.n-1
    
    t=b(k+1);
    m1=min(k,slv.band);
    for m=0:m1-1
        t=t-slv.o(k-m,m+1)*b(k-m);
    end
    b(k+1)=t;
end
for k=slv.n-1:-1:0
    t=b(k+1)/slv.d(k+1);
    m1=min(slv.n-k-1,slv.band);
    for m=0:m1-1
        t=t-slv.o(k+1,m+1)*b(k+m+2);
    end
    b(k+1)=t;
end





return
