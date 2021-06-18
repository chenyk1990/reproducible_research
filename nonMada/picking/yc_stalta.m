function [ b ] = yc_stalta( a, nsta, nlta)
% 
% nsta: short time period
% nlta: long time period
% STA/LTA
% i>=nsta,
% STA_i = (1/nsta) * sum_{j=i-nsta}^{i} a_j
% LTA_i = (1/nlta) * sum_{j=i-nlta}^{i} a_j
% i<nsta,STA_i=0,LTA_i=0
%
% Copyright: Yangkang Chen, Nov, 2016
% 
% STA/LTA = ~
% 
% Reference:
% Chen, Y., 2020, Automatic microseismic event picking via unsupervised
% machine learning, GJI, 1750–1764

sta=a;
lta=a;

sta=cumsum(sta.^2);
sta(nsta+1:end)=sta(nsta+1:end)-sta(1:end-nsta);
sta=sta/nsta;
lta=cumsum(lta.^2);
lta(nlta+1:end)=lta(nlta+1:end)-lta(1:end-nlta);
lta=lta/nlta;
sta(1:nlta-1)=zeros(nlta-1,1);
idx=find(abs(lta)<=0.0000001);
lta(idx)=lta(idx)+0.0000001;
b=sta./lta;

end


