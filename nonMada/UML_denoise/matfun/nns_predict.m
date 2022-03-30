function y = nns_predict(model,x)
% nns_predict: given NN model, x, predict y
%
%  IN   x:     features (variables)
%       model: trained model
%  OUT  y:     classfication results
%
%  Examples: pywinml/matfun/test_nn_scratch_matlab.m
%
%  Copyright (C) 2016 Yangkang Chen

[n1,n2]=size(x);
%Gradient descent parameters (I picked these by hand)
reg_lambda = 0.01;

W1=model.W1;
b1=model.b1;
W2=model.W2;
b2=model.b2;

z1=x*W1+ones(n1,1)*b1;
a1=yc_softplus(z1);
z2=a1*W2+ones(n1,1)*b2;
y=yc_softplus(z2);

return




function d1= yc_softplus(d1)
d1=log(1+exp(d1));

% d1=tanh(d1);
return
