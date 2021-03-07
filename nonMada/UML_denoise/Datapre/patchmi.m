function datami = patchmi(dataraw,patchnum,l1,l2)
% Author: Mi Zhang
% dataraw: raw data
% patchnum: the number of mini-patch
% l1: the row of mini-patch
% l2: the column of mini-patch

datami = zeros(patchnum,l1*l2);
I = dataraw;

for i = 1:patchnum
    choose_row = round((size(I,1) - l1)*rand);
    choose_col = round((size(I,2) - l2)*rand);
    I_patch = I(choose_row+1:(choose_row+l1),choose_col+1:(choose_col+l2));
    datami(i,:) = I_patch(:)';
end
end

