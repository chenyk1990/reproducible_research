function rse=rse(a,b)
%Normalized root square error (RSE) (Parekh and Selesnick, 2017)
% a and b are the vectorized estimated signal and vectorized original data, respectively. The range of RSE values is from 0 to 1. The smaller values will show better recovery quality.

rse = norm(a(:)-b(:)) / norm(b(:));
end
