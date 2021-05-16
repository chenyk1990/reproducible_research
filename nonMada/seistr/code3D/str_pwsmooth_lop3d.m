function ds=str_pwsmooth_lop3d(dn,dipi,dipx,r1,r2,eps,order)
% str_pwsmooth_lop3d: 3D plane-wave smoothing
% 
% BY Yangkang Chen, Hang Wang, and co-authors, 2019
%
% 
% INPUT:
% dn: model  noisy data
% dipi: inline slope
% dipx: xline slope
% r1,r2:       spray radius
% order:    PWD order
% eps: regularization (default:0.01);
% 
% OUTPUT:
% ds: smoothed data
% 
% Reference
% H. Wang, Y. Chen, O. Saad, W. Chen, Y. Oboue, L. Yang, S. Fomel, and Y. Chen, 2021, A Matlab code package for 2D/3D local slope estimation and structural filtering: in press.


np=(2*r1+1)*(2*r2+1);
%flattening
[u] = str_pwspray_lop3d(dn,dipi,dipx,r1,r2,order,eps);

%smoothing
ds=squeeze(sum(u,2)/np);






