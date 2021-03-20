function data = glreadsu(L,option)
fid=fopen(L,'r');
hdr=fread(fid,240,'int16');
nz=hdr(58);
fseek(fid,0,1);
pos=ftell(fid);
nx=pos/(nz*4+240);
fseek(fid,0,-1);
datap=fread(fid,nx*(nz+60),'float');
datap=reshape(datap,[nz+60 nx]);
fclose(fid);
data=datap(61:nz+60,1:nx);
if nargin>1
if(strcmp(option,'dis'))
    imagesc(data);
end
end

clear datap
end


