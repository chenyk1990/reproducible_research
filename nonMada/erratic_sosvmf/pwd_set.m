function [out] = pwd_set(adj,w,diag,offd,pp,inp)
%pwd_set: matrix multiplication
%
% INPUT:
%
% adj:adjoint flag
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% n1: trace length
% nw: PWD filter(accuracy) order (default nw=1)
% pp: slope                  (1D array)
% diag: defined diagonal .   (1D array) (useless)
% offd: defined off-diagonal (2D array) (useless)
% inp: model
% out: data 


%define PWD object(struct)
% w=struct;
% w.n=n1;
% w.na=2*nw+1;
% w.a=zeros(n1,w.na);%2D array
% w.b=zeros(w.na,1);%1D array

%apfilt_init(nw)
% nw=(w.na-1)/2;
n=w.n;
nw=(w.na-1)/2;

%% pwd_set
tmp=zeros(n,1);

if adj
    for i=0:n-1
        tmp(i+1)=0.0;
    end
    for i=0:n-1
        for j=0:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                tmp(k+1)=tmp(k+1)+w.a(k+1,j+1)*inp(i+1);
            end
            
        end
    end
    
    for i=0:n-1
        out(i+1)=0.0;
    end
    
    for i=nw:n-nw-1
        for j=0:w.na-1
            k=i+j-nw;
            out(k+1)=out(k+1)+w.a(i+1,j+1)*tmp(i+1);
        end
    end
else
    for i=0:n-1
        tmp(i+1)=0.0;
    end
    for i=nw:n-nw-1
        for j=0:w.na-1
            k=i+j-nw;
            tmp(i+1)=tmp(i+1)+w.a(i+1,j+1)*inp(k+1);
        end
    end
    for i=0:n-1
        out(i+1)=0.0;
        for j=0:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                out(i+1)=out(i+1)+w.a(k+1,j+1)*tmp(k+1);
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
