function snr=SNR(I,In)
[nt,nhx,nhy,nx,ny]=size(I);
I=reshape(I,[nt,nhx*nhy*nx*ny]);
In=reshape(In,[nt,nhx*nhy*nx*ny]);
Ps=sum(sum((I-mean(mean(I))).^2));
Pn=sum(sum((I-In).^2));
snr=10*log10(Ps/Pn);
end
