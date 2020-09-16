function d1 = yc_micro_emd(d)
% Simple EMD denoise
% by Hang Wang and Yangkang Chen
% Zhejiang University
% 
% Reference:
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51-KS61.

[n1,n2]=size(d);

d1=zeros(size(d));
for i2=1:n2
   tmp=emd(d(:,i2)); 
   d1(:,i2)=sum(tmp(2:end,:),1)'; 
end

end

