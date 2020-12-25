function [d_recon] = mssa5d_rsvd(dn0,nt,nf,nhx,nhy,nx,ny,dt,flow,fhigh,N,K,P,eps,Niter,a,mask)
%% data processing - MSSA

% Transform into F-X domain
DATA_FX=fft(dn0,nf,1);
DATA_FX0=zeros(nf,nhx,nhy,nx,ny);

MASK=squeeze(mask(1,:,:,:,:));

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1
    ihigh=floor(nf/2)+1;
end

lhx=floor(nx/2)+1;
lhxx=nx-lhx+1;
lhy=floor(ny/2)+1;
lhyy=ny-lhy+1;

lx=floor(nhx/2)+1;
lxx=nhx-lx+1;
ly=floor(nhy/2)+1;
lyy=nhy-ly+1;

% main loop on frequency
for k=ilow:ihigh

    S_obs=squeeze(DATA_FX(k,:,:,:,:)); 
    Sn_1=S_obs;
    
    for iter=1:Niter
        
        M4=P_HH(Sn_1,lx,ly,lxx,lhx,lhy,lhxx,lhyy);
        %M4=rsvd(M4,N);
        M4=rSVDbasic(M4,N,P,K);
        Sn=P_AA(M4,nf,nhx,nhy,nx,ny,lx,ly,lxx,lhx,lhy,lhxx,lhyy);
        
        Sn=a(iter)*S_obs+(1-a(iter))*MASK.*Sn+(1-MASK).*Sn;
        
        if norm(reshape(Sn,nhx*nhy,nx*ny)-reshape(Sn_1,nhx*nhy,nx*ny),'fro')<eps
            break;
        end
        
        Sn_1=Sn;
        
    end
    
    DATA_FX0(k,:,:,:,:)=reshape(Sn,1,nhx,nhy,nx,ny);
   
    if(mod(k,5)==0)
        fprintf( 'F %d is done!\n\n',k);
    end
    
end

% Back to TX (the output)
for k=nf/2+2:nf
    DATA_FX0(k,:,:,:,:) = conj(DATA_FX0(nf-k+2,:,:,:,:));
end

d_recon=real(ifft(DATA_FX0,[],1));
d_recon=d_recon(1:nt,:,:,:,:);

end