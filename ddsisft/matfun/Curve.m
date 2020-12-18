% test
clc;clear;close all;

d=zeros(251,128);
rsf_read(d,'../hyper/hyper-lap-zoom.rsf');
figure;imagesc(d);

C = fdct_wrapping_cyk(d,0,2);
% dr=ifdct_wrapping_cyk(C,0,2);dr=real(dr);
% figure;imagesc([d,dr,d-dr]);

nscale=max(size(C));
C1=C;    
    coef=[];
    for iscale=1:nscale
    	   nangle=max(size(C(iscale)));
		for iangle=1:nangle
            if iangle >floor(nangle/2);
                C1{iscale}{iangle}=zeros(size(C1{iscale}{iangle}));
            end        
       		end
    end		

ds1=ifdct_wrapping_cyk(C1,0,2);ds1=real(ds1);
figure;imagesc([d,ds1,d-ds1]);



% img = fdct_wrapping_dispcoef_cyk(C);
% figure;imagesc(abs(img));

