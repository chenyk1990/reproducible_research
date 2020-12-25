function [dout]=P_AA(din,nf,nhx,nhy,nx,ny,lx,ly,lxx,lhx,lhy,lhxx,lhyy)
% Averaging the block Hankel matrix to output the result
dout=zeros(nhx,nhy,nx,ny);
M2=P_A2(din,nhx,nhy,lx,ly,lxx,lhx,lhy,lhxx,lhyy);

for i=1:nhx
    for j=1:nhy
        dm=reshape(M2(i,j,:,:),lhx*lhy,lhxx*lhyy);
        dori=P_A(dm,nx,ny,lhx,lhy);
        dout(i,j,:,:)=dori;
    end
end
return
end

function [dout]=P_A2(din,nhx,nhy,lx,ly,lxx,lhx,lhy,lhxx,lhyy)
% Averaging the block Hankel matrix to output the result
dout=zeros(nhx,nhy,lhx*lhy,lhxx*lhyy);

for j=1:nhy
    if j<ly
        for id=1:j
            dout(:,j,:,:) =dout(:,j,:,:)+ ave_antid2(din(1+(j-1)*lx*lhx*lhy-...
                (id-1)*lx*lhx*lhy:j*lx*lhx*lhy-(id-1)*lx*lhx*lhy,1+(id-1)*lxx*lhxx*lhyy:lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy),nhx,lx,lhx,lhy,lhxx,lhyy)/j;
        end
    else
        for id=1:(nhy-j+1)           
            dout(:,j,:,:) =dout(:,j,:,:)+ ave_antid2(din((ly-1)*lx*lhx*lhy+...
                1-(id-1)*lx*lhx*lhy:ly*lx*lhx*lhy-(id-1)*lx*lhx*lhy,(j-ly)*lxx*lhxx*lhyy+1+(id-1)*lxx*lhxx*lhyy:(j-ly+1)*lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy),nhx,lx,lhx,lhy,lhxx,lhyy)/(nhy-j+1);
        end
    end
end
return
end

function [dout]=P_A(din,nhx,nhy,lx,ly)
% Averaging the block Hankel matrix to output the result
lxx=nhx-lx+1;
lyy=nhy-ly+1;
dout=zeros(nhx,nhy);

for j=1:nhy
    if j<ly
        for id=1:j
            dout(:,j) =dout(:,j)+ ave_antid(din(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx))/j;
        end
    else
        for id=1:(nhy-j+1)
            dout(:,j) =dout(:,j)+ ave_antid(din((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx))/(nhy-j+1);
        end
    end
end
return
end

function [dout] =ave_antid2(din,nhx,lx,lhx,lhy,lhxx,lhyy);
% averaging along antidiagonals
dout=zeros(nhx,lhx*lhy,lhxx*lhyy);

for i=1:nhx
    if i<lx
        for id=1:i            
            dout(i,:,:)=dout(i,:,:) + reshape(din(1+(i-id)*lhx*lhy:(i-id+1)*lhx*lhy,1+(id-1)*lhxx*lhyy:id*lhxx*lhyy)/i,1,lhx*lhy,lhxx*lhyy);
        end
    else
        for id=1:nhx+1-i
            dout(i,:,:)=dout(i,:,:) + reshape(din(1+(lx-id)*lhx*lhy:(lx-id+1)*lhx*lhy,1+(i-lx+id-1)*lhxx*lhyy:(i-lx+id)*lhxx*lhyy)/(nhx+1-i),1,lhx*lhy,lhxx*lhyy);
        end
    end
end
dout=reshape(dout,nhx,1,lhx*lhy,lhxx*lhyy);
return
end

function [dout] =ave_antid(din);
% averaging along antidiagonals
[n1,n2]=size(din);
nout=n1+n2-1;
dout=zeros(nout,1);
for i=1:nout
    if i<n1
        for id=1:i
            dout(i)=dout(i) + din(i-(id-1),1+(id-1))/i;
        end
    else
        for id=1:nout+1-i
            dout(i)=dout(i) + din(n1-(id-1),1+(i-n1)+(id-1))/(nout+1-i);
        end
    end
end
return
end

