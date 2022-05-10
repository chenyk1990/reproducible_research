function [dout]=yc_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb)
% yc_bandpass: Bandpass filtering.
%
% Aug, 05, 2020
% by Yangkang Chen
%
% INPUT
% din:      input data
% dt:       sampling
% flo:      Low frequency in band, default is 0
% fhi:      High frequency in band, default is Nyquist
% nplo=6:   number of poles for low cutoff
% nphi=6:   number of poles for high cutoff
% phase=0:  y: minimum phase, n: zero phase
% verb=0:   verbosity flag
%
% OUTPUT
% dout:     output data
%
% RFERENCE
% November 2012 program of the month: http://ahay.org/blog/2012/11/03/program-of-the-month-sfbandpass/
%
% Example
% mada/test_yc_bandpass.m
%
if size(din,1)==1 && size(din,2)>1
    din=din(:);
end

[n1,n22,n33,n44,n55]=size(din);
n2=n22*n33*n44*n55;
dout=zeros(n1,n2);

d1=dt;
eps=0.0001;

if nargin==1
    flo=0;
    fhi=0.5;
    nphi=6;
    nplo=6;
    phase=0;
    verb=1;
end

if nargin==4
    nphi=6;
    nplo=6;
    phase=0;
    verb=0;
end

if flo<0
    error('Negative flo');
else
    flo=flo*d1;
end

if fhi<0
    error('Negative flo');
else
    fhi=fhi*d1;
    if flo>fhi
        error('Need flo < fhi\n');
    end
    if 0.5<fhi
        error('Need fhi < Nyquist\n');
    end
end

if nplo<1
    nplo=1;
end
if nplo>1 && ~phase
    nplo=nplo/2;
end
if nphi<1
    nphi=1;
end
if nphi>1 && ~phase
    nphi=nphi/2;
end

if (verb)
    fprintf("flo=%g fhi=%g nplo=%d nphi=%d\n",flo,fhi,nplo,nphi);
end

if flo>eps
    blo=butter_init(0,flo,nplo);
else
    blo=[];
end
if fhi<0.5-eps
    bhi=butter_init(1,fhi,nphi);
else
    bhi=[];
end


for i2=0:n2-1
    trace=din(:,i2+1);
    
    if ~isempty(blo)
        trace=butter_apply(blo,n1,trace);
        if ~phase
            trace=reverse (n1, trace);
            trace=butter_apply (blo, n1, trace);
            trace=reverse (n1, trace);
        end
    end
    
    if ~isempty(bhi)
        trace=butter_apply(bhi,n1,trace);
        if ~phase
            trace=reverse (n1, trace);
            trace=butter_apply (bhi, n1, trace);
            trace=reverse (n1, trace);
        end
    end
    dout(:,i2+1)=trace(:);
    
end


dout=reshape(dout,n1,n22,n33,n44,n55);

return



function [bw]=butter_init(low,cutoff,nn)
% butter_init: initialize
% Aug, 5, 2020
% Yangkang Chen
% 
% INPUT
% low:      low-pass (or high-pass)
% cutoff:   cut off frequency
% nn:       number of poles
% 
% OUTPUT
% bw:       butterworth struct
% 
bw=struct;
arg=2*pi*cutoff;
sinw=sin(arg);
cosw=cos(arg);

bw.nn=nn;
bw.low=low;
bw.den=zeros(2,floor((nn+1)/2));

if mod(nn,2)>0
    if low
        fact=(1+cosw)/sinw;
        bw.den(1,floor(nn/2)+1)=1./(1.+fact);
        bw.den(2,floor(nn/2)+1)=1.-fact;
    else
        fact=sinw/(1.+cosw);
        bw.den(1,floor(nn/2)+1)=1./(fact+1.0);
        bw.den(2,floor(nn/2)+1)=fact-1.0;
    end
end

fact=das_ifnot(low,sin(0.5*arg),cos(0.5*arg));
fact=fact*fact;

for j=0:floor(nn/2)-1
    ss=sin(pi*(2*j+1)/(2*nn))*sinw;
    bw.den(1,j+1)=fact/(1.+ss);
    bw.den(2,j+1)=(1-ss)/fact;
end
bw.mid=-2.*cosw/fact;
return


function [x]=butter_apply(bw,nx,x)
% butter_apply: filter the data (in place)
% 
%Implementation is inspired by D. Hale and J.F. Claerbout, 1983, Butterworth
%dip filters: Geophysics, 48, 1033-1038.
%
% Aug, 5, 2020
% Yangkang Chen
% 
% INPUT
% bw: butterworth struct
% nx: size of x
% x: input data
% 
% OUTPUT
% x: output data
d1=bw.mid;
nn=bw.nn;

if mod(nn,2)>0
    d0=bw.den(1,floor(nn/2)+1);
    d2=bw.den(2,floor(nn/2)+1);
    x0=0;
    y1=0;
    for ix=0:nx-1
        x1=x0;x0=x(ix+1);
        y0=das_ifnot(bw.low,(x0 + x1 - d2 * y1)*d0,(x0 - x1 - d2 * y1)*d0);
        x(ix+1)=y0;
        y1=y0;
    end 
end

for j=0:floor(nn/2)-1
    d0=bw.den(1,j+1);
    d2=bw.den(2,j+1);
    x1=0;x0=0;y1=0;y2=0;
    for ix=0:nx-1
        x2=x1;x1=x0;x0=x(ix+1);
        y0=das_ifnot(bw.low,(x0 + 2*x1 + x2 - d1 * y1 - d2 * y2)*d0,(x0 - 2*x1 + x2 - d1 * y1 - d2 * y2)*d0);
        y2=y1;x(ix+1)=y0;y1=y0;
    end
end


return

function [trace]=reverse(n1,trace)
% reverse a trace (in place)
for i1=0:floor(n1/2)-1
    t=trace(i1+1);
    trace(i1+1)=trace(n1-i1);
    trace(n1-i1)=t;
end

return

