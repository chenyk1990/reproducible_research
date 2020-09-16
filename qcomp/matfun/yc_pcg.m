function [m,misfit] = yc_pcg(operator,Param,d,m0,Niter_in,Niter_out,verb)
% yc_pcg: Precondioned CG for solving sparsity-promoting inverse problems
% min || d - Fm||_2^2 + \mu || m ||_1
%
% By Yangkang Chen
% Dec, 2016
%
% INPUT
% operator: forward operator F
% Param:    parameter struct of the forward operator
% d:        RHS of the inverse problem
% m0:       initial model estimation
% Niter_in: inner iteration NO
% Niter_out:outer iteration NO
% verb:     verbosity
% 
% OUTPUT
% m:        estimated model
% misfit:   misfit history
% 
% REFERENCE
% Chen, 2018, GEO, Automatic velocity analysis using high-resolution hyperbolic Radon transform
% Wang et al., 2020, TGRS, Q-compensated denoising of seismic data

u = m0;
P = ones(size(u));
kc = 1;
Misfit = [];
m = u;
for l = 1:Niter_out
    di = feval(operator,P.*u,Param,1);
    r = d-di;
    
    g = feval(operator,r,Param,-1);
    g = g.*P;
    s = g;
    gammam = cgdot(g); %r^T_{k-1}r_{k-1}
    k = 1;
    while  k<=Niter_in
        q = feval(operator,P.*s,Param,1);
        den = cgdot(q);
        alpha = gammam/(den+1.e-8);
        u = u+alpha*s;
        r = r-alpha*q;
        misfit(kc) =  cgdot(r);
        g = feval(operator,r,Param,-1);
        g = g.*P;
        gamma = cgdot(g);%r^T_{k}r_{k}
        beta = gamma/(gammam + 1.e-7);
        gammam = gamma;
        s = g + beta*s;
        if verb
            fprintf('Iteration = %d Misfit=%0.5g\n',k,misfit(k));
        end
        k = k + 1;
        kc = kc + 1;
    end
    
    m = P.*u;
%     P=sqrt(abs(m))+0.00001; %or P=sqrt(abs(m))+0.000001; 
    % m=Pu, P=diag(abs(m0)^0.5);
    
    P = abs(m/max(abs(m(:))))+0.001; %sometimes this one is better, why?
end
return


function [out] = cgdot(in)

% Dot product

temp  =   in.*conj(in);
out = sum(temp(:));

return;

