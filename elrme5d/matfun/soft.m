function x = soft(y, T)
% x = soft(y, T)
%
% SOFT THRESHOLDING
% for real or complex data.
%
% INPUT
%   y : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%
% OUTPUT
%   x : output of soft thresholding
%
% If x and T are both multidimensional, then they must be of the same size.

x = max(1 - T./abs(y), 0) .* y;
end