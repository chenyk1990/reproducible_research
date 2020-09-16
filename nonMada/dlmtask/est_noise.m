function sigma = est_noise(DataCube)
[m,n,p]=size(DataCube);
data=reshape(DataCube,[p,m*n])';
sigma=zeros(p,1);

for i=1:p
    y=data(:,i);
    X=data;
    X(:,i)=[];
    sigma(i)=std(X*((X'*X)\(X'*y))-y);
end
