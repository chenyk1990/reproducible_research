function writdat1(name,data,dataY,dataX)
if nargin==2
    [dataY,dataX]=size(data);
end
strFileName=sprintf('%s',name);
fp=fopen(strFileName,'wb');
for i=1:dataY
    for j=1:dataX
        fwrite(fp,data(i,j),'float64');
    end
end
fclose(fp);
end


