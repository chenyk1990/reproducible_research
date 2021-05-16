function ds = str_pwsmooth_lop2d(dn,dip,ns,order,eps)
% str_pwsmooth_lop: plane-wave smoothing 
% 
% BY Yangkang Chen, Hang Wang, and co-authors, 2019
%
% INPUT:
% dn: model   noisy data
% dip: slope (2D array)
% ns:       spray radius
% order:    PWD order
% eps: regularization (default:0.01);
% ds: data   smoothed data

n1=size(dn,1);
n2=size(dn,2);

ns2=2*ns+1;%spray diameter



% w=zeros(ns2,1);
% w1=ones(n1,n2)/ns2;
% 
% for is=0:ns2-1
%     w(is+1)=1;%ns+1-abs(is-ns);
% end


% for Normalization
% t=zeros(n12,1);

utmp=str_pwspray_lop2d(dn,dip,ns,order,eps);

u=reshape(utmp,n1,ns2,n2);

% for i2=0:n2-1
%     for i1=0:n1-1
%         ws=w1(i1+1,i2+1);
%         for is=0:ns2-1
%             ds(i2*n1+i1+1)=ds(i2*n1+i1+1)+u(i1+1,is+1,i2+1)*w(is+1)*ws;
%         end
%     end
% end
% 
% ds=reshape(ds,n1,n2);
ds=squeeze(sum(u,2)/ns2);

return




