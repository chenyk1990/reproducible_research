function E_out=siga(E,delta)
[Nz,Nx]=size(E);
E_out=E;
for i=1:Nz
    for j=1:Nx
        if abs(E(i,j))<=delta
            E_out(i,j)=E(i,j);
        else
            E_out(i,j)=sign(E(i,j))*delta;
        end
    end
end
