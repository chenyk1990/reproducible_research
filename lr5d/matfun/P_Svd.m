function [dout]=P_Svd(din,N,K)
% Rank reduction on the block Hankel matrix


%      [U,D,V]=svds(din,N); % a little bit slower for small matrix
%      dout=U*D*V';
% %      
    [U,D,V]=svd(din);
%     for j=1:N
%         D(j,j)=D(j,j)*(1-D(N+1,N+1)^K/D(j,j)^K);
%     end    
    
    dout=U(:,1:N)*D(1:N,1:N)*(V(:,1:N)');

return
end