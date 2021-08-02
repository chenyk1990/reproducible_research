function [w,diag,offd] = pwd_define(forw,diag,offd,n1,nw,pp)
%pwd_define: matrix multiplication
%
% INPUT:
%
% forw:forward or backward
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% n1: trace length
% nw: PWD filter(accuracy) order (default nw=1)
% pp: slope                  (1D array)
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% ifset: only define diagonal,offdiagonal or apply pwd_set


%define PWD object(struct)
w=struct;
w.n=n1;
w.na=2*nw+1;
w.a=zeros(n1,w.na);%2D array
w.b=zeros(w.na,1);%1D array

%apfilt_init(nw)
% nw=(w.na-1)/2;
n=w.n;
nb=2*nw;

%     for (i=0; i < n; i++) {
% 	passfilter (pp[i], w->b);


%slv=sf_banded_init(n1,nb);
% diag=zeros(n1,1);
% offd=zeros(n1,nb);

for i=0:n-1
    w.b=passfilter(pp(i+1),nw);
    for j=0:w.na-1
        if(forw)
            w.a(i+1,j+1)=w.b(w.na-j);
        else
            w.a(i+1,j+1)=w.b(j+1);
        end
        
    end
end

for i=0:n-1
    for j=0:w.na-1
        k=i+j-nw;
        if k>=nw && k<n-nw
            aj=w.a(k+1,j+1);
            diag(i+1)=diag(i+1)+aj*aj;
        end
    end
    for m=0:2*nw-1
        for j=m+1:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                aj=w.a(k+1,j+1);
                am=w.a(k+1,j-m);
                offd(i+1,m+1)=offd(i+1,m+1)+am*aj;
            end
        end
    end
    
end

return


function [a,b] = passfilter(p,nw)
% passfilter: find filter coefficients£¨verfied)
% All-pass plane-wave destruction filter coefficients
%
% INPUT:
% p: slope
% a: output filter (n+1) (1D array)
%
n=nw*2;
b=zeros(n+1,1);
for k=0:n
    bk=1;
    for j=0:n-1
        if (j<n-k)
            bk=bk*(k+j+1.0)/(2*(2*j+1)*(j+1));
        else
            bk=bk*1.0/(2*(2*j+1));
        end
        
    end
    b(k+1)=bk;
end

for k=0:n
    ak=b(k+1);
    for j=0:n-1
        if j<n-k
            ak=ak*(n-j-p);
        else
            ak=ak*(p+j+1);
        end
    end
    a(k+1)=ak;
end
return
