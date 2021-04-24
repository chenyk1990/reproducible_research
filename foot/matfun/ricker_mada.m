function [ dout ] = ricker_mada(din,dt,freq,deriv)
% Frequency-domain filtering
%
% Yangkang Chen, Zhejiang University, June, 2019
%
% INPUT:
%
% dt: time interval
% freq: peak frequency for Ricker wavelet (in Hz)
% deriv: apply a half-order derivative filter (1: apply; 2: not apply)
% din: input data (column-wise) (nt*ntrace)
%
% OUTPUT:
% dout:output data
%
%  ALSO SEE
%  ~/chenyk/matlibcyk/mada/test_ricker_mada.m
%
%  Reference:
%  http://ahay.org/blog/2012/12/23/program-of-the-month-sfhalfint/ (expression for integration is wrong!!!)


%determine frequency sampling (for real to complex FFT)

%setting parameters for ricker_init();

if nargin==3
    deriv=0;
end
[n1,n22,n33,n44,n55]=size(din);

[n1,n2]=size(din);
nfft = 2^nextpow2(n1);
freq=freq*dt;

if deriv==1
    order=2;
else
    order=0;
end

nw=nfft/2+1;
dw=1./(nfft*freq);
shape=zeros(nw,1);

for iw=0:nw-1
    w=iw*dw;
    w=w*w;
    
    switch order
        
        case 2 %half-order derivative
            cw=sqrt(2*pi/nfft+i*iw*2*pi);
            shape(iw+1)=cw*w*exp(1-w);
        case 0
            shape(iw+1)=w*exp(1-w);
        otherwise
            error('Wrong order');
            %   matlab fft is different from Mada version (so no need to dividing cf by nfft)
            
    end
    
    
end


for i2=1:n2
    
    dout(:,i2)=freqfilt(din(:,i2),nfft,shape);
    %or using the operator form
    % [ ~, dout(:,i2) ] = freqfilt_lop(0,0,shape, n1,n1, din(:,i2), [] );
end


if n2~=n22 && n2==n22*n33*n44*n55
    dout=reshape(dout,n1,n22,n33,n44,n55);
end
return


function [ dout ] = freqfilt(din,n,shape)

n1=size(din,1);

nw=n/2+1;
tmp=fft(din,n,1);
tmp(1:nw)=tmp(1:nw).*shape;

for iw=nw+1:n
    tmp(iw)=conj(tmp(n-iw+2));
end

tmp2=real(ifft(tmp,[],1));%inverse fft
dout=tmp2(1:n1);

return




