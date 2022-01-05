function [DATA_f] = fxarma(DATA,flow,fhigh,dt,lf,mu,verb);
%FXARMA: SNR enhancement using fx projection filtering (FX ARMA).
%
%  [DATA_f] = fxarma(DATA,dt,lf,mu,flow,fhigh);
% 
%  IN   DATA:   the data matrix, columns are traces
%       dt:     sampling interval in sec
%       lf:     lenght of operator (lenght of the filter)
%       mu:     pre-whitening in %
%       flow:   min  freq. in the data in Hz
%       fhigh:  max  freq. in the data in Hz
% 
%  OUT  DATA_f: filtered data 
%
%  Reference: Sacchi and Kuehl, 2000, FX ARMA filters, 70.th. Ann. Internat. 
%             Mtg., Soc. Expl. Geophys., Expanded Abstracts
%       
%  Example: see testfxarma.m
%
%  Copyright (C) 2014, Texas Consortium for Computational Seismology
%  Author: Yangkang Chen, The University of Texas at Austin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%


 [nt,ntraces] = size(DATA);
 nf = 2^nextpow2(nt);
 
 DATA_FX_f = zeros(nf,ntraces);
 DATA_FX_b = zeros(nf,ntraces);

% First and last samples of the DFT.

 ilow  = floor(flow*dt*nf)+1; 

  if ilow<1; 
   ilow=1; 
  end;

 ihigh = floor(fhigh*dt*nf)+1;

  if ihigh > floor(nf/2)+1; 
   ihigh=floor(nf/2)+1; 
  end

% Transform to FX

 DATA_FX = fft(DATA,nf,1);

 for k = ilow:ihigh;
  aux_in  = DATA_FX(k,:)';
  aux_out = arma_modeling(aux_in,lf,mu);
  DATA_FX(k,:) = aux_out';
    if(mod(k,5)==0 && verb==1)
        fprintf( 'F %d is done!\n\n',k);
    end
 end;

% Honor symmetries

 for k=nf/2+2:nf
  DATA_FX(k,:) = conj(DATA_FX(nf-k+2,:));
 end

% Back to TX (the output) 

 DATA_f = real(ifft(DATA_FX,[],1));
 DATA_f = DATA_f(1:nt,:);

return

function [y] = arma_modeling(x,lf,mu);
%ARMA_MODELING: autoregressive moving averaging (ARMA) modeling of 1D spatial data
%
%  IN    x:   data 
%        lf:  length of the operator
%        mu:  pre-whitening parameter
%      
%  OUT   yf:  prediction of the data using AR modeling
% 
%  Copyright (C) 2014, Yangkang Chen
%
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%

   nx = length(x);
   y=zeros(nx,1);
   C=convmtx(x,lf);     % convolution matrix
   R=C'*C/nx;
   [g,Pw]=eigs(R,1,'SM'); 
   g0=g(1);
   g=g/g0;              % filter (ARMA)
   e=C*g;   
   G=convmtx(g,nx);
   
   
   B = G'*G;  beta = B(1,1)*mu; 
   w = (B + beta*eye(nx))\G'*e;    % Estimated additive noise
   y=x-w;

return
