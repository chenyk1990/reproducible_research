function [dout]=P_Lmafit(din,N)
% Rank reduction on the block Hankel matrix
    [m,n]=size(din);
    [data,position]=PreServ(din);    
    opts=[];
    [X,Y,~] = lmafit_mc_adp(m,n,N,position,data,opts);
    dout=X*Y;
return
end

function [data, position]=PreServ(input)
[m n]=size(input);
% count=1;
% 
% for i=1:n
%     
%     for j=1:m
%         
%         if input(j,i)~=0
%             
%             data(count)=input(j,i);            
%             position(count)=(i-1)*m+j;
%             count=count+1;
%             
%         end
%         
%     end
% 
% end

position=1:m*n;
data=reshape(input,1,m*n);

end