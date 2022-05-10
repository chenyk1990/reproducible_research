function [dip] = yc_dip2dmask(din,mask,niter,liter,order,eps_dv, eps_cg, tol_cg,rect,verb)
% yc_dip2d_i: dip estimation based on shaping regularized PWD algorithm (verified)
% (independent implementation)
% 
% exactly the same with dip2d_shan_GN.M
% 
% BY Yangkang Chen, Nov, 05, 2019
% INPUT
% din: input data (nt*nx)
% niter: number of nonlinear iterations
% liter: number of linear iterations (in divn)
% order: accuracy order
% eps_dv: eps for divn  (default: 0.01)
% eps_cg: eps for CG    (default: 1)
% tol_cg: tolerence for CG (default: 0.000001)
% rect:  smoothing radius (ndim*1)
% verb: verbosity flag
%
% OUTPUT
% dip:  2D slope
% 
% Example:
% test/test_PWD_yc_dip_i.m
if nargin==1
    mask=[];
    niter = 5;
    liter = 20;
    eps_dv = 0.01;
    eps_cg = 1;
    tol_cg = 0.000001;
%     tol_cg = 0.0;
    order=2;
    rect(1) = 5;
    rect(2) = 5;
    rect(3) = 1;
    verb=1;%debug
end

if nargin==2
    niter = 5;
    liter = 20;
    eps_dv = 0.01;
    eps_cg = 1;
    tol_cg = 0.000001;
%     tol_cg = 0.0;
    order=2;
    rect(1) = 5;
    rect(2) = 5;
    rect(3) = 1;
    verb=1;%debug
end

if size(din) ~= size(mask)
    mask=[];
    niter=mask;
    liter=niter;
    order=liter;
    eps_dv=order;
    eps_cg=eps_dv;
    tol_cg=eps_cg;
    rect=tol_cg;
    verb=rect;
end

p0 = 0.0;
nj1 = 1;
nj2 = 1;

dim = 3;

n = zeros(dim, 1);
n1 = size(din,1);
n2 = size(din,2);
n(1) = n1;
n(2) = n2;
n(3) = 1;

n123 = n(1) * n(2) * n(3);

ratio=zeros(size(din));
dip=zeros(size(din));
mask0=mask;
for i1=1:n(1)
    for i2=1:n(2)-1
            mask(i1,i2)= (mask(i1,i2) && din(i1,i2+1));
    end
end
% figure;imagesc([mask0,mask,mask0-mask]);colorbar;
% figure;imagesc([mask0-mask]);colorbar;

for iter=1:niter
[u1,u2] = conv_allpass(din,dip,order); % linearization using the updated dip

if ~isempty(mask)
%     u1(find(mask==0))=0;
%     u2(find(mask==0))=0;
    u1=u1.*mask;
    u2=u2.*mask;
end
[ ratio ] = yc_divne(-u2, u1, liter, rect, n, eps_dv, eps_cg, tol_cg,verb);
dip=dip+ratio;
end

return

function [u1,u2] = conv_allpass(din,dip,order)
% Convolutional operator implemented by an allpass filter
%
% Linearized inverse problem
% C'(\sigma)d\Delta \sigma =  C(\sigma)d
%
% OUTPUT:
% u1: C'(\sigma)d (denominator)
% u2: C(\sigma)d  (numerator)

u1=zeros(size(din));
u2=zeros(size(din));

[n1,n2]=size(din);

if order==1
    nw=1;
else
    nw=2;
end
filt1=zeros(2*nw+1,1);
filt2=zeros(2*nw+1,1);

for i1=nw+1:n1-nw
    for i2=1:n2-1
        
        if order==1
            filt1=B3d(dip(i1,i2));
            filt2=B3(dip(i1,i2));
        else
            filt1=B5d(dip(i1,i2));
            filt2=B5(dip(i1,i2));
        end
        for iw=-nw:nw
            u1(i1,i2)=u1(i1,i2)+(din(i1+iw,i2+1)-din(i1-iw,i2))*filt1(iw+nw+1);
            u2(i1,i2)=u2(i1,i2)+(din(i1+iw,i2+1)-din(i1-iw,i2))*filt2(iw+nw+1);
        end
    end
end

return

function [b3 ] = B3(sigma)
% B3 coefficient
% sigma: slope

b3(1)=(1-sigma)*(2-sigma)/12;
b3(2)=(2+sigma)*(2-sigma)/6;
b3(3)=(1+sigma)*(2+sigma)/12;

return

function [b3d ] = B3d(sigma)
% B3 coefficient derivative
% sigma: slope
b3d(1)=-(2-sigma)/12-(1-sigma)/12;
b3d(2)=(2-sigma)/6-(2+sigma)/6;
b3d(3)=(2+sigma)/12+(1+sigma)/12;

return

function [b5 ] = B5(sigma)
% B5 coefficient
% sigma: slope

b5(1)=(1-sigma)*(2-sigma)*(3-sigma)*(4-sigma)/1680;
b5(2)=(4-sigma)*(2-sigma)*(3-sigma)*(4+sigma)/420;
b5(3)=(4-sigma)*(3-sigma)*(3+sigma)*(4+sigma)/280;
b5(4)=(4-sigma)*(2+sigma)*(3+sigma)*(4+sigma)/420;
b5(5)=(1+sigma)*(2+sigma)*(3+sigma)*(4+sigma)/1680;

return

function [b5d ] = B5d(sigma)
% B5 coefficient derivative
% sigma: slope

b5d(1)=-(2-sigma)*(3-sigma)*(4-sigma)/1680-...
    (1-sigma)*(3-sigma)*(4-sigma)/1680-...
    (1-sigma)*(2-sigma)*(4-sigma)/1680-...
    (1-sigma)*(2-sigma)*(3-sigma)/1680;

b5d(2)=-(2-sigma)*(3-sigma)*(4+sigma)/420-...
    (4-sigma)*(3-sigma)*(4+sigma)/420-...
    (4-sigma)*(2-sigma)*(4+sigma)/420+...
    (4-sigma)*(2-sigma)*(3-sigma)/420;

b5d(3)=-(3-sigma)*(3+sigma)*(4+sigma)/280-...
    (4-sigma)*(3+sigma)*(4+sigma)/280+...
    (4-sigma)*(3-sigma)*(4+sigma)/280+...
    (4-sigma)*(3-sigma)*(3+sigma)/280;

b5d(4)=-(2+sigma)*(3+sigma)*(4+sigma)/420+...
    (4-sigma)*(3+sigma)*(4+sigma)/420+...
    (4-sigma)*(2+sigma)*(4+sigma)/420+...
    (4-sigma)*(2+sigma)*(3+sigma)/420;

b5d(5)=(2+sigma)*(3+sigma)*(4+sigma)/1680+...
    (1+sigma)*(3+sigma)*(4+sigma)/1680+...
    (1+sigma)*(2+sigma)*(4+sigma)/1680+...
    (1+sigma)*(2+sigma)*(3+sigma)/1680;

return

