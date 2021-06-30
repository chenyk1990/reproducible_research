function A=h_hankel(fa,Nx,Lx)
%Construct the hankel matrix
%A=hankel(fa,50,40), A=hankel(fa,51,38)
if nargin==2
    Lx=1;
end
[n,m]=size(fa);
Ny=m+1-Nx;
Ly=n+1-Lx;
for i=1:n
    for j=1:Ny     
        H(i,(j-1)*Nx+1:j*Nx)=fa(i,j:Nx+j-1);%Hankel matrix
    end
end
for k=0:(Ly-1)
    for i=1:Lx
        for j=1:Ny
            A((j+k*Ny),(i-1)*Nx+1:i*Nx)=H(i+k,(j-1)*Nx+1:j*Nx);
        end
    end
end
end