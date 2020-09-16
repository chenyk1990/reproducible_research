function [A, S] = MTSNMF(X, r, lambda, opts)
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
%                                achieved when C(t)-C(t-1) < C(t)*opts.tol.
%                                Default is 1e-5.
%               opts.verbose - Verbose or not. Default is true.
%
% Outputs:
%    A      - K-by-1 cell array containing K M-by-R basis matrices.
%    S      - R-by-N coefficient matrix shared by all tasks.


error(nargchk(3, 4, nargin, 'struct'));

maxit = 3000;
tol = 1e-5;
alpha = 1;
tau = 0.2;
verbose = true;

if exist('opts', 'var')
    if isfield(opts, 'maxit')
        maxit = opts.maxit;
    end
    if isfield(opts, 'tol')
        tol = opts.tol;
    end
    if isfield(opts, 'tau')
        tau = opts.tau;
    end
    if isfield(opts, 'alpha')
        alpha = opts.alpha;
    end
    if isfield(opts, 'verbose')
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


if isa(X{1}, 'double')
    eps_val = 1e-16;
elseif isa(X{1}, 'single')
    eps_val = single(1e-8);
else
    error('Invalid input data type.');
end



K = length(X);
A = cell(size(X));
[n, m] = size(X{1});
S = max(rand(r, m, class_type), eps_val);

for k = 1:K
    A{k} = max(rand(n, r, class_type), eps_val);
end

rho_S = ((2*n+1)*(m+r)*K + 2*(r+1)*m)/(2*m*(r+4));
rho_A = (2*n*m+2*m*r+2*n*r+n)/(n*(2*r+7));
% rho_A = 10;
lambda_sum = sum(cell2mat(lambda));

S_row_max = max(S, [], 2);
[max_S, max_S_row_idx] = max(S_row_max);


A_col_max = zeros(K, r, class_type);

for k=1:K
    A_col_max(k, :) = max(A{k}, [], 1);
end

[max_A, max_A_col_idx] =  max(A_col_max, [], 2);

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
        idx = randperm(r);
        for j_count = 1:r
            j = idx(j_count);
            W = S(j,:)+(U(j,:)-V(j,:)*S-lambda_sum)/V(j,j);
%             W = S(j,:)+(U(j,:)-V(j,:)*S)/V(j,j);
            W = max(W, eps_val*max_S);
            S(j,:) = W;
            S_row_max(j) = max(W);
            if S_row_max(j) >= max_S
                max_S = S_row_max(j);
                max_S_row_idx = j;
            elseif j == max_S_row_idx
                [max_S, max_S_row_idx] = max(S_row_max);
            end
        end
%         S = wthresh(S,'h',lambda_sum/(10000));
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
            idx = randperm(r);
            for j_count = 1:r
                j = idx(j_count);
                W = A{k}(:,j)+(U(:,j)-A{k}*V(:,j))/V(j,j);
                W = max(W, eps_val*max_A(k));
                A{k}(:,j) = W;
                A_col_max(k, j) = max(W);
                if A_col_max(k, j) >= max_A(k)
                    max_A(k) = A_col_max(k, j);
                    max_A_col_idx(k) = j;
                elseif j == max_A_col_idx(k)
                    [max_A(k), max_A_col_idx(k)] = max(A_col_max(k, :));
                end
            end
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
    C = 1/2*C+lambda_sum*sum(sum(S));
    
    if verbose
        fprintf('Iteration %d.\tObjective function value C=%f\n', itrcount, C);
        toc(tStart);
    end
    
    if abs(C_old-C) < tol*C
        break;
    end
end

end

