function dout = phase_correction_cyk(din,c);
%PHASE_CORRECTION: Apply a constant phase correction to seismic data.
%              
%  [dout] = phase_correction(din,c);
%
%  IN   din:      Seismic data (traces in columns);
%       c:        Constant Phase rotation in degrees
%
%  OUT  dout:     data after correction
%
%
%  Reference: Longbottom, J., Walden, A.T. and White, R.E. (1988) Principles and 
%             application of maximum kurtosis phase estimation. Geophysical 
%             Prospecting, 36, 115-138.
%
%  Copyright (C) 2008,  Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%          Yangkang Chen, Zhejiang University
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
%  Test example:
%  test/test_phase_rotate_ricker.m
% 

  c = c*pi/180;
  [nt,nx]=size(din);
  nf = 2^nextpow2(nt);
  [nt,nx]=size(din);

  Din = fft(din,nf,1);

  Phase = exp(-i*[0;-c*ones(nf/2-1,1);0;c*ones(nf/2-1,1)]);

   for k=1:nx;
     Din(:,k) = Din(:,k).*Phase;
   end

 dout = ifft(Din,nf,1);
 dout = real(dout(1:nt,:));

 return;
