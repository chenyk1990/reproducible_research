function [dout]=MA(D, e, ws)
%The moving-average (MA) filter

b = (1/ws)*ones(1,ws);
dout = filter(b,e,D);
end
