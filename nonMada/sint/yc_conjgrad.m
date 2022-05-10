function [ x ] = yc_conjgrad(opP,opL,opS, p, x, dat, eps_cg, tol_cg, N,ifhasp0,par_P,par_L,par_S,verb)
% yc_conjgrad: conjugate gradient with shaping
%
% BY Yangkang Chen, Nov, 04, 2019
% Modified by Yangkang Chen, Nov, 09, 2019 (fix the "adding" for each oper)
%
% INPUT
% opP: preconditioning operator
% opL: forward/linear operator
% opS: shaping operator
% p: preconditioned model
% x: estimated model
% dat:  data
% N:  number of iterations
% eps_cg:  scaling
% tol_cg:  tolerance
% ifhasp0: flag indicating if has initial model
% par_P: parameters for P
% par_L: parameters for L
% par_S: parameters for S
% verb: verbosity flag
%
% OUPUT
% m: estimated model
%

np=length(p(:));
nx=par_L.nm;    %model size
nd=par_L.nd;    %data size

if ~isempty(opP)
    d=-dat; %nd*1
    r=feval(opP,d,par_P,0,0);
else
    r=-dat;
end



if ifhasp0
    x=feval(opS,p,par_S,0,0);

    if ~isempty(opP)
        d=feval(opL,x,par_L,0,0);
        par_P.d=r;%initialize data
        r=feval(opP,d,par_P,0,1);
    else
        par_L.d=r;%initialize data
        r=feval(opL,x,par_L,0,1);
    end
else
    p=zeros(np,1);%define np!
    x=zeros(nx,1);%define nx!
end

dg=0;
g0=0;
gnp=0;
r0=sum(r.*r);   %nr*1

for n=1:N
    gp=eps_cg*p; %np*1
    gx=-eps_cg*x; %nx*1

    if ~isempty(opP)
        d=feval(opP,r,par_P,1,0);%adjoint
        par_L.m=gx;%initialize model
        gx=feval(opL,d,par_L,1,1);%adjoint,adding
    else
        par_L.m=gx;%initialize model
        gx=feval(opL,r,par_L,1,1);%adjoint,adding
    end

    par_S.m=gp;%initialize model
    gp=feval(opS,gx,par_S,1,1);%adjoint,adding
    gx=feval(opS,gp,par_S,0,0);%forward,adding
    if ~isempty(opP)
        d=feval(opL,gx,par_P,0,0);%forward
        gr=feval(opP,d,par_L,0,0);%forward
    else
        gr=feval(opL,gx,par_L,0,0);%forward
    end
    
    gn = sum(gp.*gp); %np*1
    
    if n==1
        g0=gn;
        sp=gp; %np*1
        sx=gx; %nx*1
        sr=gr; %nr*1
    else
        alpha=gn/gnp;
        dg=gn/g0;
        
        if alpha < tol_cg || dg < tol_cg
            if verb
                fprintf("convergence in %d iterations, alpha=%g, gd=%g\n",n,alpha,dg);
            end
            break;
        end
        
        gp=alpha*sp+gp;
        t=sp;sp=gp;gp=t;
        
        gx=alpha*sx+gx;
        t=sx;sx=gx;gx=t;
        
        gr=alpha*sr+gr;
        t=sr;sr=gr;gr=t;
    end
    
    beta=sum(sr.*sr)+eps_cg*(sum(sp.*sp)-sum(sx.*sx));    
    if verb
        fprintf('iteration: %d res: %f grad: %f!\n',n,sum(r.* r) / r0,dg);
    end
    
    alpha=-gn/beta;
    
    p=alpha*sp+p;
    x=alpha*sx+x;
    r=alpha*sr+r;
    
    gnp=gn;
end


return