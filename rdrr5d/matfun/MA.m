function [dout]=MA(D, e, ws)
%The moving-average (MA) filter
%
%  KEY REFERENCE
%  Oboue et al., 2021, Robust damped rank-reduction method for simultaneous denoising and reconstruction of 5-D seismic data, Geophysics, 86, V71–V89.

b = (1/ws)*ones(1,ws);
dout = filter(b,e,D);
end
