function  z = SoftTh(s,thld)
% SoftTh : regularized iterative soft thresholding operator
%      s : Signal
%   thld : regularized iterative soft thresholding parameter
%      z : Signal approximated

        z = sign(s).*max(0,abs(s)- abs(thld));
end
