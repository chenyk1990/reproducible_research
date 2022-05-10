function [w1] = pwsmooth_set(dip,n1,n2,ns,order,eps)
%pwsmooth_set 

ns2=2*ns+1;%spray diameter
n12=n1*n2;

u=zeros(n1,ns2,n2);

w=zeros(ns2,1);
w1=zeros(n1,n2);

for is=0:ns2-1
    w(is+1)=ns+1-abs(is-ns);
end

% for Normalization
t=zeros(n12,1);


for i1=0:n12-1
    w1(i1+1)=1.0;
end

w1tmp=w1(:);
[w1tmp,t] = pwsmooth_lop(0,0,dip,w1,n1,n2,ns,order,eps,n12,n12,w1tmp,t);

w1=reshape(w1tmp,n1,n2);

for i1=0:n12-1
    if 0~= t(i1+1)
        w1(i1+1)=1.0/t(i1+1);
    else
        w1(i1+1)=0.0;
    end
    
end


end

