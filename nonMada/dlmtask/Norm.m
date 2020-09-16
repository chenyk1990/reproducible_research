function [NDnorm,min_val,max_val] = Norm(ND)

min_val=min(min(min(ND)));
max_val=max(max(max(ND)));
NDnorm=(ND-min_val)/(max_val-min_val);

return