function b = banded_solve(n,band,diag,offd,b)
%banded_solve: Banded matrix solver
%
% Yangkang Chen, Zhejiang University, Oct, 2018
%
% INPUT:
%
% n:    matrix size
% band: band size
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
%
% OUTPUT:
% b: solution
%
% Test:
%  A=[1,0.5,0.1,0;0.5,2,0.3,0.05;0.1,0.3,3.0,0.2;0,0.05,0.2,4.0];
%  b=[1,2,1.5,1.0]';
%  x1=A\b;
%  diag=[1,2,3,4];offd=[[0.5;0.3;0.2;0],[0.1;0.05;0;0]];
%  x2=banded_solve(4,2,diag,offd,b);
%  norm(x1-x2)
% 
% Reference:
% Fomel, Madagascar, api/c/banded.c

%define Band object(struct)
slv=struct;
slv.n=n;
slv.band=band;
slv.d=zeros(n,1);       %1D array
slv.o=zeros(n-1,band);  %2D array


% define the banded matrix
for k=0:slv.n-1
    t=diag(k+1);
    m1=min(k,slv.band);
    for m=0:m1-1
        t=t-slv.o(k-m,m+1)*slv.o(k-m,m+1)*slv.d(k-m);
    end
    slv.d(k+1)=t;
    n1=min(slv.n-k-1,slv.band);
    for n=0:n1-1
        t=offd(k+1,n+1);
        m1=min(k,slv.band-n-1);
        for m=0:m1-1
            t=t-slv.o(k-m,m+1)*slv.o(k-m,n+m+2)*slv.d(k-m);
        end
        slv.o(k+1,n+1)=t/slv.d(k+1);
    end
    
    
end

% the solver 
for k=1:slv.n-1
    
    t=b(k+1);
    m1=min(k,slv.band);
    for m=0:m1-1
        t=t-slv.o(k-m,m+1)*b(k-m);
    end
    b(k+1)=t;
end
for k=slv.n-1:-1:0
    t=b(k+1)/slv.d(k+1);
    m1=min(slv.n-k-1,slv.band);
    for m=0:m1-1
        t=t-slv.o(k+1,m+1)*b(k+m+2);
    end
    b(k+1)=t;
end



end