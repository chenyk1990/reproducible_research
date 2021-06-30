function data_out=hwl_double_MSSA_noise(data_in,flow,fhigh,dt,k,a,A_iter,lambda)
%Dealiased and denoised MSSA
%dt: time interval
%k: rank
%a: interpolation interval
%A_iter: iteration numbers

[n,m,l]=size(data_in);
data_out=zeros(n,m,l);
f_data_in=data_in;
Nx=floor(m/2)+1;
Lx=floor(l/2)+1;
Ny=m+1-Nx;
Ly=l+1-Lx;
Nx2=m*(a-1)+Nx;
Lx2=l*(a-1)+Lx;
%Ly2=l+1-Lx;
num1=m;
%fa_re=zeros(l,m);
t=n;x=m;z=l;

nf=2^nextpow2(t);

% Transform into F-X domain
f_data_in=fft(data_in,nf,1);
%fshotc=fft(shotc,nf,1);
data=zeros(t,x*a,z*a);

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1;
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1;
    ihigh=floor(nf/2)+1;
end

G=zeros(m*a,m);
GG=zeros(m*a,m);
for i=1:m
    G(1+(i-1)*a,i)=1;
    GG(2+(i-1)*a:i*a,i)=ones(a-1,1);
end

G2=zeros(l*a,l);
GG2=zeros(l*a,l);
for i=1:l
    G2(1+(i-1)*a,i)=1;
    GG2(2+(i-1)*a:i*a,i)=ones(a-1,1);
end

for i=ilow:ihigh
    fa_re=zeros(l,m);
    i_low=round(i/a);
    temp_low=f_data_in(i_low,:,:);%这里fa是频率切片
    fa_low=reshape(temp_low,[m,l]);
    fa_low=fa_low';
    A_low=h_hankel(fa_low,Nx,Lx);
    [U_low,S_low,V_low]=svd(A_low);
    
    
    temp=f_data_in(i,:,:);
    fa=reshape(temp,[m,l]);
    fa_new1=G*fa;
    fa_new1=fa_new1';
    fa_new=G2*fa_new1;
    fa_iter=fa_new;
    for A_i=1:A_iter
        A=h_hankel(fa_iter,Nx2,Lx2);
        M=U_low(:,1:k)*U_low(:,1:k)'*A;%TSVD
        fa=anti_Hankel_2(M,Nx2,Ny,Lx2,Ly,m*a,l*a);
        fa_re=fa;
        for alj=1:a:l*a
            for ali=1:a:m*a
                fa_re(alj,ali)=0;
                %   fa_re(:,ali)=fa(:,1+(ali-1)*a);
            end
        end
        % figure;imagesc(abs(fa_re))
        if A_i==1
            ra=0;
        else
            ra=(1/(A_iter-A_i+1)).^lambda;
        end
        fa_iter=(ra)*(fa-fa_re)+(1-ra)*fa_new+fa_re;
    end
    
    data(i,:,:)=reshape(fa_iter',[1,m*a,l*a]);
    fprintf('f=%d is done\n',i);
end


for k=nf/2+2:nf
    data(k,:,:) = conj(data(nf-k+2,:,:));
end
data_out=real(ifft(data,[],1));
data_out=data_out(1:t,:,:);





