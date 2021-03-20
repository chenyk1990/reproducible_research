function   glwritesu( L,Dt,path,optional)
[Nz Nx]=size(L);

zero1=zeros(1,57);
zero2=zeros(1,61);
if nargin==1
    Dt=0.001;
    Dt=Dt*1000*1000;
    fid=fopen('data.su','w+');
    for i=1:Nx
         fwrite(fid,zero1,'int16');
         fwrite(fid,Nz,'int16');
         fwrite(fid,Dt,'int16');
         fwrite(fid,zero2,'int16');
         fwrite(fid,L(:,i),'float'); 
    end
    fclose(fid);
end

if nargin==2
    Dt=Dt*1000*1000;
    fid=fopen('data.su','w+');
    for i=1:Nx
         fwrite(fid,zero1,'int16');
         fwrite(fid,Nz,'int16');
         fwrite(fid,Dt,'int16');
         fwrite(fid,zero2,'int16');
         fwrite(fid,L(:,i),'float'); 

    end
    fclose(fid);
end

if nargin==3
    Dt=Dt*1000*1000;
    fid=fopen(path,'w+');
   
  
    for i=1:Nx
         fwrite(fid,zero1,'int16');
         fwrite(fid,Nz,'int16');
         fwrite(fid,Dt,'int16');
         fwrite(fid,zero2,'int16');
         fwrite(fid,L(:,i),'float'); 
    end
    fclose(fid);
end

if nargin==4
    Dt=Dt*1000*1000;
    fid=fopen(path,'a+');
    for i=1:Nx
         fwrite(fid,zero1,'int16');
         fwrite(fid,Nz,'int16');
         fwrite(fid,Dt,'int16');
         fwrite(fid,zero2,'int16');
         fwrite(fid,L(:,i),'float'); 

    end
    fclose(fid);
end


end

