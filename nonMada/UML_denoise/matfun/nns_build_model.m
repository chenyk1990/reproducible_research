function model = nns_build_model(X,y,nn_hdim,niter,lrate,nepoch,verb,param)
%  nns_build_model: training a NN
%  IN   X:  training dataset (N*m)
%       y:  given classification mode (N*1)
%       nn_hdim: hidden layer of the NN (L) -> W1:m*L;z1=N*L;a1=N*L;W2:L*c;z2=N*c;a2(y_pre)=N*c;
%       niter: training iterations (Niter)
%       verb: debugging verbosity
%
%  OUT  model: trained model
%
%  Examples: pywinml/matfun/test_nn_scratch_matlab.m
%
%  Copyright (C) 2016 Yangkang Chen

[n1,n2]=size(X);
nn_input_dim=size(X,2);

% Gradient descent parameters (I picked these by hand)
if nargin==7
    epsilon = 0.01; % learning rate for gradient descent
    reg_lambda = 0.01; % regularization strength
    param=struct;
else
    epsilon = param.eps;
    reg_lambda=param.lambda;
end

if isfield(param, 'nc') %number of classes
    nn_output_dim=param.nc;
else
    nn_output_dim=size(X,2);
end

if isfield(param, 'vf') %verbosity frequency (every 10,100 or 1000 iterations)
    vf=param.vf;
else
    vf=100;
end

randn('state',201617);
W1=randn(nn_input_dim,nn_hdim)/sqrt(nn_input_dim);
randn('state',201617);
W2=randn(nn_hdim,nn_output_dim)/sqrt(nn_output_dim);
% max(W1(:))
b1=zeros(1,nn_hdim);
b2=zeros(1,nn_output_dim);

% size(W1) %W1: [64,128]
% size(W2) %W2: [128,64]
% size(b1) %b1: [1,128]
% size(b2) %b2: [1,64]
% size(X)  %X:  [7047,64]
model=struct('W1',W1,'b1',b1,'W2',W2,'b2',b2);

rand('state',2020);
inds=randperm(size(X,1));
[~,inds2]=sort(inds,'ascend');
XX=X(inds,:);
yy=y(inds,:);
% X=XX(inds2,:);%reconstruct
% y=yy(inds2,:);
% nepoch=70;
patch=floor(size(X,1)/nepoch);%batch size
for ie=1:nepoch
    if ie*patch<=7047
        X=XX(1+(ie-1)*patch:ie*patch,:);
        y=yy(1+(ie-1)*patch:ie*patch,:);
    else
        X=XX(1+(ie-1)*patch:end,:);
        y=yy(1+(ie-1)*patch:end,:);
    end
    n1=size(X,1);
    for iter=1:niter
        
        % forward propagation
        %      size(X)
        %      size(W1)
        %      size(b1)
        z1=X*W1+ones(n1,1)*b1;%size: z1 [7047,128]
        a1=yc_softplus(z1);   %size: a1 [7047,128]
        %     a1=z1;
        z2=a1*W2+ones(n1,1)*b2;%size:z2 [7047,64]
        yp=yc_softplus(z2);   %size: yp [7047,64]
        %     yp=z2;
        %      size(W2)
        %      size(b2)
        
        % backpropagation
        %     delta3=probs;
        %     for i1=1:n1
        %     delta3(i1,y(i1)+1)=delta3(i1,y(i1)+1)-1;
        %     end
        %     size(delta3)
        %     size(y)
        %     size(probs)
        %     size(y)
        delta3=(yp-y).*yc_sigmoid(z2); %delta3: [7047,64]
        %     delta3
        %     max(delta3(:))
        dW2=a1'*delta3;    %a1': [128,7047], dW2: [128,64]
        db2=sum(delta3,1); %[1,64]
        %     size(db2)
        delta2=delta3*W2'.*yc_sigmoid(z1);%[7047,64] x [64,128] .* [7047,128] = [7047,128]
        dW1=X'*delta2;  %[64,7047] x [7047,128] = [64,128]
        db1=sum(delta2,1);%[1,128]
        
        % add regularization terms (b1 and b2 does not require regularization)
        dW2=dW2+reg_lambda*W2;
        dW1=dW1+reg_lambda*W1;
        
        %gradient descent parameter update
        %     epsilon0=epsilon*((niter-iter)/(niter-1)).^0.3;
        epsilon0=epsilon;
        epsilon0=lrate;
        W1=W1-epsilon0*dW1;
        b1=b1-epsilon0*db1;
        W2=W2-epsilon0*dW2;
        b2=b2-epsilon0*db2;
        if 1
            %         epsilon0
            %     max(yp(:))
            %     max(dW1(:))
            %     [max(dW1(:)),max(dW2(:)),max(db1(:)),max(db2(:))]
        end
        
        %     [max(dW1(:)),max(dW2(:))]
        %     max(a1(:))
        %     max(delta3(:))
        %     max(z2(:))
        
        model.W1=W1;
        model.b1=b1;
        model.W2=W2;
        model.b2=b2;
        
%         if verb==1 && mod(iter,vf)==0
%             %      if verb==1
%             fprintf('Loss after epoch %d/%d, iteration %d/%d: %f\n',ie,nepoch,iter, niter, nns_calculate_loss(X,y,model,reg_lambda));
%         end
    end
        if verb==1 
            %      if verb==1
            fprintf('Loss after epoch %d/%d is %f \n',ie,nepoch, nns_calculate_loss(X,y,model,reg_lambda));
        end
end

return


function d1= yc_softplus(d1)
d1=log(1+exp(d1));

% d1=tanh(d1);
return

function d1= yc_sigmoid(d1)
d1=1./(1+exp(-d1));
% d1=1-tanh(d1).^2;
return
