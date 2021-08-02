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









