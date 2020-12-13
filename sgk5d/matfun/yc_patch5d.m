function [ X ] = yc_patch5d( A,mode,l1,l2,l3,l4,l5,s1,s2,s3,s4,s5)
%decompose 4D/5D data into patches:
%
% by Yangkang Chen
% March, 2020
%
% Input:
%   D: input image
%   mode: patching mode
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
%             April 2, 2020 (3D-5D)
%
% Examples:
%    ~/test/test_yc_ksvd_denoise3d.m
%    ~/test/test_yc_ksvd_denoise5d.m
%    ~/test/test_yc_ksvd_recon5d.m
% 
% References:
% Wang et al., 2020, Fast dictionary learning for high-dimensional seismic reconstruction, IEEE Transactions on Geoscience and Remote Sensing, doi: 10.1109/TGRS.2020.3030740
% Zhou et al., 2020, Statistics-guided dictionary learning for automatic coherent noise suppression, IEEE Transactions on Geoscience and Remote Sensing, in press
% Saad, O. and Y. Chen, 2020, PATCHUNET: A fully-unsupervised and highly-generalized deep learning approach for random noise suppression, Geophysical Prospecting. 
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

[n1,n2,n3,n4,n5]=size(A);

if mode==1 %possible for other patching options
    
    if nargin==2
        l1=4;l2=4;l3=4;l4=4;l5=4;s1=2;s2=2;s3=2;s4=2;s5=2;
    end
    
    if nargin==7
        s1=round(l1/2);s2=round(l2/2);s3=round(l3/2);s4=round(l4/2);s5=round(l5/2);
    end
    
    tmp=mod(n1-l1,s1);
    if tmp~=0
        A=[A;zeros(s1-tmp,n2,n3,n4,n5)];
    end
    
    tmp=mod(n2-l2,s2);
    if tmp~=0
        A=[A,zeros(size(A,1),s2-tmp,n3,n4,n5)];
    end
    
    tmp=mod(n3-l3,s3);
    if tmp~=0
        A=cat(3,A,zeros(size(A,1),size(A,2),s3-tmp,n4,n5));%concatenate along the third dimension
    end
    
    tmp=mod(n4-l4,s4);
    if tmp~=0
        A=cat(4,A,zeros(size(A,1),size(A,2),size(A,3),s4-tmp,n5));%concatenate along the forth dimension
    end    

    tmp=mod(n5-l5,s5);
    if tmp~=0
        A=cat(5,A,zeros(size(A,1),size(A,2),size(A,3),size(A,4),s5-tmp));%concatenate along the fifth dimension
    end    
    
    [N1,N2,N3,N4,N5]=size(A);
    X=[];
    for i1=1:s1:N1-l1+1
        for i2=1:s2:N2-l2+1
            for i3=1:s3:N3-l3+1
                for i4=1:s4:N4-l4+1
                    for i5=1:s5:N5-l5+1
                        tmp=reshape(A(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1,i4:i4+l4-1,i5:i5+l5-1),l1*l2*l3*l4*l5,1);
                        X=[X,tmp];
                    end
                end
            end
        end
    end
    
end




end

