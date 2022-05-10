function [dout] = yc_mask_lop(din,par,adj,add )
% yc_mask_lop: Mask operator
% BY Yangkang Chen, Aug, 10, 2021
% INPUT
% d: model/data
% par: parameter (mask,nm,nd)
% adj: adj flag
% add: add flag
% OUTPUT
% m: data/model

mask=par.mask;
nm=length(mask(:));
nd=nm;
par.nm=nm;
par.nd=nd;
if adj==1
    d=din;
    if isfield(par,'m') && add==1
        m=par.m;
    else
        m=zeros(par.nm,1);
    end
else
    m=din;
    if isfield(par,'d') && add==1
        d=par.d;
    else
        d=zeros(par.nd,1);
    end
end

[ m,d ] = yc_adjnull( adj,add,nm,nd,m,d );

for im=1:nm
    if mask(im)
        if adj==1
            m(im)=m(im)+d(im);
        else %forward
            d(im)=d(im)+m(im);
        end
    end
end

if adj==1
    dout=m;
else
    dout=d;
end

return

