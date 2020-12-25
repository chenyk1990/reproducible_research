function [dout]=P_HH(din,lx,ly,lxx,lhx,lhy,lhxx,lhyy)

H2=P_H2(din,lhx,lhy,lhxx,lhyy);
dout=P_H4(H2,lx,ly,lxx,lhx,lhy,lhxx,lhyy);

return
end

function [dout]=P_H4(din,lx,ly,lxx,lhx,lhy,lhxx,lhyy)
% forming block Hankel matrix
[~,nhy,~,~]=size(din);

for j=1:nhy
    h2=squeeze(din(:,j,:,:));
    r=P_H3(h2,lhx,lhy,lhxx,lhyy);
    if j<ly
        for id=1:j
            dout(1+(j-1)*lx*lhx*lhy-(id-1)*lx*lhx*lhy:j*lx*lhx*lhy-(id-1)*lx*lhx*lhy,1+(id-1)*lxx*lhxx*lhyy:lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy) = r;
        end
    else
        for id=1:(nhy-j+1)
            dout((ly-1)*lx*lhx*lhy+1-(id-1)*lx*lhx*lhy:ly*lx*lhx*lhy-(id-1)*lx*lhx*lhy,(j-ly)*lxx*lhxx*lhyy+1+(id-1)*lxx*lhxx*lhyy:(j-ly+1)*lxx*lhxx*lhyy+(id-1)*lxx*lhxx*lhyy)=r;
        end
    end
end
return
end

function [dout]=P_H3(din,lhx,lhy,lhxx,lhyy)
% forming 3rd order block Hankel matrix

[nhx,~,~]=size(din);
lx=floor(nhx/2)+1;

for j=1:nhx
    if j<lx
        for id=1:j
            dout(1+(j-1)*lhx*lhy-(id-1)*lhx*lhy:j*lhx*lhy-(id-1)*lhx*lhy,1+...
                (id-1)*lhxx*lhyy:lhxx*lhyy+(id-1)*lhxx*lhyy) = din(j,:,:);
        end
    else
        for id=1:(nhx-j+1)
            dout((lx-1)*lhx*lhy+1-(id-1)*lhx*lhy:lx*lhx*lhy-(id-1)*lhx*lhy,...
                (j-lx)*lhxx*lhyy+1+(id-1)*lhxx*lhyy:(j-lx+1)*lhxx*lhyy+(id-1)*lhxx*lhyy)=din(j,:,:);
        end
    end
end

return
end

function [dout]=P_H2(din,lhx,lhy,lhxx,lhyy)
% forming block Hankel matrix
[nhx,nhy,~,~]=size(din);

dout=zeros(nhx,nhy,lhx*lhy,lhxx*lhyy);

for ii=1:nhx
    for jj=1:nhy
        d=squeeze(din(ii,jj,:,:));
        dout(ii,jj,:,:)=P_H(d,lhx,lhy);
    end
end
return
end

function [dout]=P_H(din,lx,ly)
% forming block Hankel matrix
[nhx,nhy]=size(din);
lxx=nhx-lx+1;
lyy=nhy-ly+1;

for j=1:nhy
    r=hankel(din(1:lx,j),[din(lx:nhx,j)]);
    if j<ly
        for id=1:j
            dout(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx) = r;
        end
    else
        for id=1:(nhy-j+1)
            dout((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx)=r;
        end
    end
end
return
end


