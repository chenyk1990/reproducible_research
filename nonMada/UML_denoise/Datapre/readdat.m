function datanew=readdat(name,dataY,dataX)
datanew=zeros(dataX,dataY);
strFileName=sprintf('%s',name);
fid = fopen(strFileName, 'r');
for i=1:dataX
    for j=1:dataY
        datanew(i,j)=fread(fid,1,'float');
    end
end
datanew=datanew';
fclose(fid);


