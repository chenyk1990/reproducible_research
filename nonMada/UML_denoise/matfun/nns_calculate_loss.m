function loss = nns_calculate_loss(X,y,model,lambda)
%  nns_calculate_loss: calculate loss function
%  IN   X:  samples
%       y:  given classification mode
%       model: trained NN model
%  OUT  loss: loss function value
%
%  Examples: pywinml/matfun/test_nn_scratch_matlab.m
%
%  Copyright (C) 2016 Yangkang Chen

[n1,n2]=size(X);

%Gradient descent parameters (I picked these by hand)
reg_lambda = lambda;

W1=model.W1;
b1=model.b1;
W2=model.W2;
b2=model.b2;

z1=X*W1+ones(n1,1)*b1;
a1=yc_softplus(z1);
z2=a1*W2+ones(n1,1)*b2;
yp=yc_softplus(z2);

loss= sum((yp(:)-y(:)).^2);
return

function d1= yc_softplus(d1)
d1=log(1+exp(d1));
% d1=tanh(d1);
return

