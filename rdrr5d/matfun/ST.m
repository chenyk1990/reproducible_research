function [dout]=ST(Z, K)
%The soft thresholding (ST) operator 
%
%  REFERENCE
%  Oboue et al., 2021, Robust damped rank-reduction method for simultaneous denoising and reconstruction of 5-D seismic data, Geophysics, 86, V71–V89.

dout = SoftTh(Z,K); 
end
