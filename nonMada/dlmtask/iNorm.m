function [DCOut] = iNorm(DC,min_val,max_val)

DCOut=max(DC,0); DCOut=min(DC,1);
DCOut=DCOut*(max_val-min_val)+min_val;

return

