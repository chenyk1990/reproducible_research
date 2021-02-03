function checkDimensions(m,n,x,mode)

if mode == 0, return; end; % OK


if (size(x,1) == 1) && (size(x,2) == 1)
  if (m ~= 1) || (n ~= 1)
     error('Operator-scalar multiplication not yet supported');
  end
end


if mode == 1
  if (size(x,1) ~= n)
    error('Incompatible dimensions');
  end;
  
  if (size(x,2) ~= 1)
    error('Operator-matrix multiplication not yet supported');
  end;
else
  if (size(x,1) ~= m)
    error('Incompatible dimensions');
  end;
  
  if (size(x,2) ~= 1)
    error('Operator-matrix multiplication not yet supported');
  end;
end
