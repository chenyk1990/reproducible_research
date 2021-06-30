function fa=anti_Hankel_2(M,Nx,Ny,Lx,Ly,n,m)
a=Nx-Ny;
b=Lx-Ly;
M_new=zeros(Ny*Ly,Nx*Lx+(Ly-1)*Nx);
temp=zeros(Ny,Nx,m);
JS=ones(Ly,Lx);
JS_new=zeros(Ly,Lx+Ly-1);
for i=1:Ly
    M_new(1+(i-1)*Ny:i*Ny,1+(i-1)*Nx:Nx*Lx+(i-1)*Nx)=M(1+(i-1)*Ny:i*Ny,:);
    JS_new(i,i:Lx+i-1)=JS(i,:);
end
JS_1=sum(JS_new);
for i=1:Lx+Ly-1
    sump=zeros(Ny,Nx);
    for j=1:Ly
    sump=sump+M_new(1+(j-1)*Ny:j*Ny,1+(i-1)*Nx:i*Nx);
    end
    temp(:,:,i)=sump/JS_1(1,i);
end

N_new=zeros(Ny,Nx+Ny-1);
JS=ones(Ny,Nx);
JS_new=zeros(Ny,Nx+Ny-1);
for i=1:m
    w=temp(:,:,i);
    for j=1:Ny
      JS_new(j,j:Nx+j-1)=JS(j,:); 
      N_new(j,j:Nx+j-1)=w(j,:);
    end
    fa(i,:)=sum(N_new)./sum(JS_new);
end

end