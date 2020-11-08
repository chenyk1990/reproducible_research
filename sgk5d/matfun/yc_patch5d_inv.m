function [ A ] = yc_patch5d_inv( X,mode,n1,n2,n3,n4,n5,l1,l2,l3,l4,l5,s1,s2,s3,s4,s5)
% insert patches into the 4D/5D data
%
% by Yangkang Chen
% April, 2020
%
% Input:
%   D: input image
%   mode: patching mode
%   n1: first dimension size
%   n1: second dimension size
%   n3: third dimension size
%   n4: forth dimension size
%   n5: fifth dimension size
%   l1: first patch size
%   l2: second patch size
%   l3: third patch size
%   l4: fourth patch size
%   l5: fifth patch size    (when n5=1, l5=1, s5=0)
%   s1: first shifting size
%   s2: second shifting size
%   s3: third shifting size
%   s4: fourth shifting size
%   s5: fifth shifting size (when n5=1, l5=1, s5=0)
%   
% Output:
%   X: patches
%
% Modified on Dec 12, 2018 (the edge issue, arbitrary size for the matrix)
% 		      Dec 31, 2018 (tmp1=mod(n1,l1) -> tmp1=mod(n1-l1,s1))
%             Marich, 31, 2020, 2D->3D
%             April 2, 2020 (3D-5D)
%
% Examples:
%    ~/test/test_yc_ksvd_denoise3d.m
%    ~/test/test_yc_ksvd_denoise5d.m
%    ~/test/test_yc_ksvd_recon5d.m
% 
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, doi: 10.1109/TGRS.2020.3030740
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., 2020, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 222, 1717-1727.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2020, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 222, 1846?1863. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51-KS61.
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.

%% patch size l1*l2*l3*l4*l5
%l1=4;l2=4;l3=4;l4=4;l5=4;
%

if mode==1 %possible for other patching options
    
    if nargin==7
        l1=4;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;
    end
    
    if nargin==12
        s1=round(l1/2);s2=round(l2/2);s3=round(l3/2);s4=round(l4/2);s5=round(l5/2);
    end
    
    tmp1=mod(n1-l1,s1);
    tmp2=mod(n2-l2,s2);
    tmp3=mod(n3-l3,s3);
    tmp4=mod(n4-l4,s4);
    tmp5=mod(n5-l5,s5);
    
    if tmp1~=0 && tmp2~=0 && tmp3~=0 && tmp4~=0 && tmp5~=0
        A=zeros(n1+s1-tmp1,n2+s2-tmp2,n3+s3-tmp3,n4+s4-tmp4,n5+s5-tmp5);
        mask=zeros(n1+s1-tmp1,n2+s2-tmp2,n3+s3-tmp3,n4+s4-tmp4,n5+s5-tmp5);
    end
    
    if tmp1~=0 && tmp2==0 && tmp3==0 && tmp4==0 && tmp5==0
        A=zeros(n1+s1-tmp1,n2,n3,n4,n5);
        mask=zeros(n1+s1-tmp1,n2,n3,n4,n5);
    end
    
    if tmp1==0 && tmp2~=0 && tmp3==0 && tmp4==0 && tmp5==0
        A=zeros(n1,n2+s2-tmp2,n3,n4,n5);
        mask=zeros(n1,n2+s2-tmp2,n3,n4,n5);
    end
    
    if tmp1==0 && tmp2==0 && tmp3~=0 && tmp4==0 && tmp5==0
        A=zeros(n1,n2,n3+s3-tmp3,n4,n5);
        mask=zeros(n1,n2,n3+s3-tmp3,n4,n5);
    end

    if tmp1==0 && tmp2==0 && tmp3==0 && tmp4~=0 && tmp5==0
        A=zeros(n1,n2,n3,n4+s4-tmp4,n5);
        mask=zeros(n1,n2,n3,n4+s4-tmp4,n5);
    end
    
    if tmp1==0 && tmp2==0 && tmp3==0 && tmp4==0 && tmp5~=0
        A=zeros(n1,n2,n3,n4,n5+s5-tmp5);
        mask=zeros(n1,n2,n3,n4,n5+s5-tmp5);
    end
    
    if tmp1==0 && tmp2==0  && tmp3==0 && tmp4==0 && tmp5==0
        A=zeros(n1,n2,n3,n4,n5);
        mask=zeros(n1,n2,n3,n4,n5);
    end
    [tmp1,tmp2,tmp3,tmp4,tmp5]
    [N1,N2,N3,N4,N5]=size(A);
    id=0;
    for i1=1:s1:N1-l1+1
        for i2=1:s2:N2-l2+1
            for i3=1:s3:N3-l3+1
                for i4=1:s4:N4-l4+1
                    for i5=1:s5:N5-l5+1
                        id=id+1;
                        A(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1,i4:i4+l4-1,i5:i5+l5-1)=A(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1,i4:i4+l4-1,i5:i5+l5-1)+reshape(X(:,id),l1,l2,l3,l4,l5);
                        mask(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1,i4:i4+l4-1,i5:i5+l5-1)=mask(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1,i4:i4+l4-1,i5:i5+l5-1)+ones(l1,l2,l3,l4,l5);
                    end
                end
            end
        end
    end
    
    A=A./mask;
    
    A=A(1:n1,1:n2,1:n3,1:n4,1:n5);
end

return

