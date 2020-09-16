function [A, S] = MTSNMF_MU(X, r, LA, lambda, opts)
%MTSNMF multi-task sparse NMF.
%
% Description:
%    Multi-task sparse NMF solves the problem:
%    (A_1,...,A_K,S) = argmin sum_k (||X_k-A_k*S||_F^2 + lambda_k*||S||_1).
%
%    [A, S] = MTSNMF(X, r, lambda) performs a multi-task sparse NMF with
%    rank R on all M-by-N data matrices in the K-by-1 cell array X,
%    resulting in K M-by-R basis matrices included in K-by-1 cell array A
%    and a common R-by-N coefficient matrix S.
%
% Syntax:
%    [A, S] = MTSNMF(X, r, lambda, opts)
%
% Inputs:
%    X      - K-by-1 cell array containing K M-by-N data matrices to be
%               factorized.
%    r      - Rank of approximation (number of basis atoms in A_1,...,A_K)
%    lambda - Regularization parameter. K-by-1 vector containing lambda_1,
%               lambda_2,...,lambda_K.
%    opts   - Options:
%               opts.maxit   - Maximal iteration number. Default is 3000.
%               opts.tol     - Convergence tolerance. The convergence is
%                                achieved when F(t)-F(t-1) < F(t)*opts.tol.
%                                Default is 1e-5.
%               opts.verbose - Verbose or not. Default is true.
%
% Outputs:
%    A      - K-by-1 cell array containing K M-by-R basis matrices.
%    S      - R-by-N coefficient matrix shared by all tasks.
%

error(nargchk(4, 5, nargin, 'struct'));

maxit = 3000;
tol = 1e-5;
alpha = 1;
tau = 0.2;
verbose = true;

if exist('opts', 'var')
    if isfield(opts, 'maxit')
        maxit = opts.maxit;
    end
    if  isfield(opts, 'tol')
        tol = opts.tol;
    end
    if  isfield(opts, 'tau')
        tau = opts.tau;
    end
    if  isfield(opts, 'alpha')
        alpha = opts.alpha;
    end
    if  isfield(opts, 'verbose')
        verbose = opts.verbose;
    end
end



if ~iscell(lambda)
    if numel(lambda) == 1
        temp = lambda;
        lambda = cell(size(X));
        lambda(:) = {temp};
    else
        lambda = num2cell(lambda);
    end
elseif numel(lambda) == 1
    temp = lambda{1};
    lambda = cell(size(X));
    lambda(:) = {temp};
end
if numel(lambda) ~= length(lambda) || length(lambda) ~= length(X)
    error('Size of lambda is invalid.');
end

if verbose
    tStart = tic;
end

class_type = class(X{1});

K = length(X);
A = cell(size(X));
[n, m] = size(X{1});
S = rand(r, m, class_type);
S = max(S, eps(S));

for k = 1:K
    A{k} = rand(n, r, class_type);
    A{k} = max(A{k}, eps(A{k}));
end
% 
rho_A = (n*m+m*r+n*r+n)/(n*r+4*n);
rho_S = ((2*n+1)*(m+r)*K+(2*r+3)*m)/((2*r+9)*m);
% rho_A = 30;
% rho_S = 60;

sum_lambda = sum(cell2mat(lambda));

C = inf;

for itrcount = 1:maxit
    U = zeros(r, m, class_type); %sum(A'X)
    V = zeros(r, r, class_type); %sum(A'A)
    for k = 1:K
        U = U+A{k}'*X{k};
        V = V+A{k}'*A{k};
    end
    
    S_outer_old = S;
    
    for inner_s_itrcount = 1:floor(1+alpha*rho_S)
        S_inner_old = S;
        S = S.*U./(V*S+sum_lambda+eps(U));
        if norm(S-S_inner_old, 'fro') < tau*norm(S-S_outer_old, 'fro')
            break;
        end
    end
    
    V = S*S';
    for k = 1:K
        A_outer_old = A{k};
        U = X{k}*S';
        for inner_a_itrcount = 1:floor(1+alpha*rho_A)
            A_inner_old = A{k};
            A{k} = A{k}.*U./(A{k}*V+LA*A{k}+eps(U));
            if norm(A{k}-A_inner_old, 'fro') < tau*norm(A{k}-A_outer_old, 'fro')
                break;
            end
        end
    end
    
    C_old = C;
    C = 0;
    for k = 1:K
        C = C+norm(X{k}-A{k}*S, 'fro')^2;
    end
    C = 1/2*C+sum_lambda*sum(sum(S));
    
    if verbose
        fprintf('Iteration %d.\tObjective function value C=%f\n', itrcount, C);
        toc(tStart);
    end
    if C_old-C < tol*C
        break;
    end
end

end

