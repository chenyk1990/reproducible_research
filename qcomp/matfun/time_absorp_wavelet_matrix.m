function [FW,TW] = time_absorp_wavelet_matrix_real(nt,dt,fdom,Q);

%  IN   nt:     lenght of the trace
%       dt:       sampling interval in secs
%       fdom:     dominant frequency in Hz 
%       Q:        quality factor
%
%  OUT  TW:  time_absorp_wavelet_matrix

sign = mod(nt,2);
if(sign == 1)
    df = 1/(nt*dt);  
else
   nt=nt+1;
   df = 1/(nt*dt);   
end 
 
% create a time vector
  t=-ceil((nt-1)/2)*dt:dt:floor((nt-1)/2)*dt;
  
% create the Ricker wavelet
  pf=pi^2*fdom^2;
  w=(1-2.*pf*t.^2).*exp(-pf*t.^2);
  w=w./max(abs(w));
  
% the wavelet in frequency domain
  nn=ceil((nt-1)/2)+1;
  ww=[w(nn:end),w(1:nn-1)];
  W=fft(ww);
  W(1,1)=0;
  
  f_ref=10*fdom;
  
% the wavelet absorption in frequency domain
  A=zeros(nn,nt);
  for i=2:nn
      for k=1:nt
          A(i,k)=exp(-(pi*df*(i-1)*(k-1)*dt)/Q(k))...
                     *exp(j*2*log((i-1)*df/f_ref)*(i-1)*df*(k-1)*dt/Q(k));
      end
  end 
  
   FW=A;
    
  A1=A(2:nn,:);   
  A1=conj(A1);
  A1=flipud(A1);
  WA=[A;A1]; 
      
% translate to the time domain
  TWA=zeros(nt,nt);
  for i=1:nt
      TWA(:,i)=ifft(WA(:,i));
  end
      
  TW=zeros(nt,nt);
  ss=[TWA(nn+1:end,:);TWA(1:nn,:)];
  for i=1:nn
      TW(1:nn+i-1,i)=ss(nn-i+1:end,i);
  end
  for i=1:nn-1
      TW(i+1:end,i+nn)=ss(1:end-i,i+nn);
  end 
  
  if(sign == 1)
      return;
  else
      TW=TW(1:end-1,1:end-1);
  end
  

return;