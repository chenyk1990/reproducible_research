function mask = genmask(u, r, type,seed)
%GENMASK:Generate Random Sampling Mask
%
%	 mask = genmask(u,r,type)
%    u, image
%    r, data KNOWN ratio
%    type: data lose type
%   'r': random lose rows
%   'c': random lose columns
%   'p': random lose pixel
%   'seed': seed of random number generator
%
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details:
%  http://www.gnu.org/licenses/


[m,n] = size(u);
mask = zeros(m,n);

switch type
    case 'r'
        row = rperm(m,seed);
        k = fix(r*m);
        row = row(1:k);
        mask(row,:) = 1;
    case 'c'
        column = rperm(n,seed);
        k = fix(r*n);
        column = column(1:k);
        mask(:, column) = 1;
    case 'p'
        pix = rperm(m*n,seed);
        r = fix(r*m*n);
        pix = pix(1:r);
        mask(pix) = 1;
end

