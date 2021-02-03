function [Shat,relmse_hat,mse_hat] = yc_optshrink_damp(Xtil,r, K)
% Damped optshrink by Yangkang Chen, Oct, 2017
% 
% K: damping factor
% 
% Input: Xtil = Signal-plus-noise matrix
%        r    = Estimate of (effective) rank of signal matrix
% Note: r can be estimated from eyeballing the "knee" of the singular value plot
%       or from your favorite random matrix theory based test
%
% Output: Shat = Rank r denoised estimate of signal matrix
%         relmse_hat = Estimate of normalized MSE
%         mse_hat    = Estimate of MSE 
%
% Idea: Shat has same singular vectors as Xtil but "shrunk" singular values
%       The shrinkage is computed in an optimal manner from the "noise portion" of the
%       singular value spectrum
%       
% ValueOfInformation type Metric: relmse_hat quantifies 'noisiness'; if  close to 1 then 
%                                 Shat is likely very unreliable
%
% Reference: Algorithm 1 in http://arxiv.org/abs/1306.6042
% Author: Raj Rao Nadakuditi (rajnrao@umich.edu)
% Date: July 22, 2013
%
% Nested sub-functions: estimateDz and estimateDpz
% 
% Example:
%
% n = 100; m = 150; theta= 2;
% u = sign(randn(n,1))/sqrt(n); v = sign(randn(m,1))/sqrt(m);
% S = theta*u*v'; 
% Xtil = S + randn(n,m)/sqrt(m);
% [Shat,relmse_hat,mse_hat] = optshrink(Xtil,1);
% [Util,Stil,Vtil] = svds(Xtil,1);
% Stsvd = Stil*Util*Vtil';
% disp('New algorithm error vs Truncated SVD error')
% [norm(S-Shat,'fro')^2 norm(S-Stsvd,'fro')^2]
% Improvement = (1-ans(1)/ans(2))*100;
% ['Improvement = ' num2str(Improvement) ' %']


[Util,Stil,Vtil] = svd(Xtil); stil = diag(Stil);

sigmahats = stil(1:r);
Xnoise_est = Stil(r+1:end,r+1:end);

for idx = 1 : r,
    theta_hats(idx) = sqrt(1/estimateDz(Xnoise_est,sigmahats(idx)));
    Dpz(idx) = estimateDpz(Xnoise_est,sigmahats(idx));
    wopt_hats(idx) = -2/(theta_hats(idx)^2*Dpz(idx))*(1-(stil(r+1)/stil(idx))^K);
end

Shat = Util(:,1:r)*diag(wopt_hats,0)*Vtil(:,1:r)';

relmse_hat = 1 - sum(wopt_hats.^2)/sum(theta_hats.^2);
mse_hat = relmse_hat*sum(theta_hats.^2);

end

function Dz = estimateDz(X,z)

% Estimates the D-transform as in Eq. 15-(a) of http://arxiv.org/abs/1306.6042
% Assume w.l.o.g X is n x m and diagonal and that z > max(svd(X))

[n,m] = size(X);

In = eye(n); Im = eye(m);

z2IXXt = z^2*In - X*X';
z2IXtX = z^2*Im - X'*X;
invz2XtX = inv(z2IXtX);
invz2XXt = inv(z2IXXt);

D1z = 1/n*trace(z*invz2XXt);
D2z = 1/m*trace(z*invz2XtX);

Dz = D1z*D2z;

%% Equivalently Dz = (1/n*trace(z*inv(z^2*In - X*X')))*(1/m*trace(z*inv(z^2*Im-X'*X)));

end

function Dpz = estimateDpz(X,z)

% Estimates the D-transform derivative as in Eq. 15-(b) of http://arxiv.org/abs/1306.6042
% Assume w.l.o.g X is n x m and diagonal and that z > max(svd(X))

[n,m] = size(X);
In = eye(n); Im = eye(m);

z2IXXt = z^2*In - X*X';
z2IXtX = z^2*Im - X'*X;
invz2XtX = inv(z2IXtX);
invz2XXt = inv(z2IXXt);

D1z = 1/n*trace(z*invz2XXt);
D2z = 1/m*trace(z*invz2XtX);


D1zp = 1/n*trace(-2*z^2*invz2XXt^2+invz2XXt);
D2zp = 1/m*trace(-2*z^2*invz2XtX^2+invz2XtX);

Dpz = D1z*D2zp+D1zp*D2z;

end



