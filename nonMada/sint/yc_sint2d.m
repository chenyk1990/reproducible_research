function [dout] = yc_sint2d(din,mask,dip,niter,ns,order,eps)
% yc_sint2: sparse data interpolation in 2-D using shaping regularization and structural filtering.
%
% By Yangkang Chen
% Aug 10, 2021
%
% INPUT
% din: input data
% mask: sampling matrix (as always in interpolation problem)
% dip:  slope
% niter: number of iterations
% ns:    smoothing radius
% order: accuracy order
% eps:  regularization
%
% OUTPUT
% dout: interpolated data
% 
% References
% [1] Chen, Y., Chen, X., Wang, Y. and Zu, S., 2019. The interpolation of sparse geophysical data. Surveys in Geophysics, 40(1), pp.73-105.
% [2] Wang, H., Chen, Y., Saad, O.M., Chen, W., Oboué, Y.A.S.I., Yang, L., Fomel, S. and Chen, Y., 2022. A Matlab code package for 2D/3D local slope estimation and structural filtering. Geophysics, 87(3), pp.F1–F14.
% [3] Huang, G., Chen, X., Li, J., Saad, O.M., Fomel, S., Luo, C., Wang, H. and Chen, Y., 2021. The slope-attribute-regularized high-resolution prestack seismic inversion. Surveys in Geophysics, 42(3), pp.625-671.
% [4] Huang, G., Chen, X., Luo, C. and Chen, Y., 2020. Geological structure-guided initial model building for prestack AVO/AVA inversion. IEEE Transactions on Geoscience and Remote Sensing, 59(2), pp.1784-1793.


if nargin==5
    order=1;
    eps=0.01;
end

[n1,n2,n3]=size(din);
n12=n1*n2;

for i3=1:n3
    mm=din(:,:,i3);
    
    %figure out scaling and make known data mask
    lam=0;
    for ii=0:n12-1
        if mask(ii+1)~=0
            lam=lam+1;
        end
    end
    lam=sqrt(lam/n12);

    %for mask operator
    par_L=struct;
    par_L.mask=mask;
    par_L.nm=n12;
    par_L.nd=n12;
    
    %for pwsmooth operator
    [w1] = pwsmooth_set(dip,n1,n2,ns,order,eps);
    par_S=struct;
    par_S.dip=dip;
    par_S.w1=w1;
    par_S.ns=ns;
    par_S.order=order;
    par_S.eps=eps;
    par_S.nm=n12;
    par_S.nd=n12;

    eps_cg=lam*lam;
    tol_cg=0.0000019209;
    ifhasp0=1;
    verb=1;
    [mm ] = yc_conjgrad([], @yc_mask_lop, @yc_pwsmooth_lop, mm(:), mm(:), mm(:), eps_cg, tol_cg, niter,ifhasp0,[],par_L,par_S,verb);
    dout(:,:,i3)=reshape(mm,n1,n2);

end

return







