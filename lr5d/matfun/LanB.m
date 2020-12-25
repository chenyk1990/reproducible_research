function [T_new]=LanB(T,k)

[m,n]=size(T);

randn('state',201314);
z1=randn(n,1);

Q=zeros(n,k);
U=zeros(m,k);
B=zeros(k,k);

alpha=zeros(k,1);
beta=zeros(k-1,1);

Q(:,1)=z1/norm(z1,2);
y1=T*Q(:,1);
alpha(1)=norm(y1,2);
U(:,1)=y1/alpha(1);

for i=1:k-1
    
    zi1=T'*U(:,i)-alpha(i)*Q(:,i);
    for j=i:-1:1
        zi1=zi1-(zi1'*Q(:,j))*Q(:,j);
    end
    
    beta(i)=norm(zi1,2);
    Q(:,i+1)=zi1/beta(i);
    yi1=T*Q(:,i+1)-beta(i)*U(:,i);
    for j=i:-1:1
        yi1=yi1-(yi1'*U(:,j))*U(:,j);
    end
    alpha(i+1)=norm(yi1,2);
    U(:,i+1)=yi1/alpha(i+1);
    
end

for i=1:k
    B(i,i)=alpha(i);
    if(i<k)
        B(i,i+1)=beta(i);
    end
end

T_new=U*B*Q';

return