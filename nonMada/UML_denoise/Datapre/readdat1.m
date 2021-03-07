function datanew=readdat1(name,dataY,dataX)
datanew=zeros(dataY,dataX);
strFileName=sprintf('%s',name);
fid = fopen(strFileName, 'r');
for i=1:dataY
    for j=1:dataX
        datanew(i,j)=fread(fid,1,'float64');
    end
end
fclose(fid);


