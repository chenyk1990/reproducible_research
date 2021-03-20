% DEMO script of erratic noise suppression using iterative curvelet thresholding
% 
% References
% Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
% Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
% Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805?1823.


clc;
clear;

is_real=1;           % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
finest=2;            % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
alpha=1.2;           % 噪音标准差的alpha倍阈值（1.2左右较为理想）

%% %输入测试模型
% infile='~/chenyk.data2/various/cyksmall/zq_rc_input.su';
infile='zq_rc_input.su';
[ori]=glreadsu(infile,'su');
ori=ori(300:2:end,:);
Input=ori;
[Nt,Nx]=size(Input);
%% % 加载奇异单道噪音以及高斯随机噪音
mask=rand(1,Nx);
mask(logical(mask<0.9))=0;
mask(logical(mask>=0.9))=1;
for i=1:Nt
    Input(i,:)=Input(i,:)+5*10^-9*randn(1,Nx).*mask;
end
Input=Input+4*10^-10*randn(Nt,Nx);
glwritesu(Input,0.001,'ori.su')
%% %曲波阈值滤波函数
F=ones(Nt,Nx);                                  % ones(n)返回n*n的1矩阵，频率域
X=fftshift(ifft2(F))*sqrt(prod(size(F)));       % prod返回size(F)的乘积,X是一个见脉冲，标准差为1
C=fdct_wrapping(X,0,finest);                    % 进行曲波变换用的是复曲波,最外层为小波变换
% Compute norm of curvelets (exact)
E=cell(size(C));
for s=1:length(C)
    E{s}=cell(size(C{s}));
    for w=1:length(C{s})
        A=C{s}{w};       
        E{s}{w}=sqrt(sum(sum(A.*conj(A)))/prod(size(A)));    %计算A的模，计算标准差，复数域单点相乘
    end
end
%% %确定阈值
CInput=fdct_wrapping(Input,is_real,finest);     %进行曲波变换用的是实曲波,最外层为小波变换
Smax=length(CInput);
Sigma=alpha*median(median(abs(CInput{Smax}{1})))/0.58;     %求取噪音标准差，选取最大尺度
sigma=[Sigma,5*Sigma,2*Sigma, Sigma, 0.6*Sigma,Sigma/5];
Sigma=sigma(1);     
%% %初次去噪处理，得到初始去噪结果（该步得到的便是常规利用曲波基去噪的结果）
Ct=CInput;
for s=2:length(CInput)
    thresh=Sigma+Sigma*s;    %最外层设置为4*sigma
    for w=1:length(CInput{s})
        Ct{s}{w}=CInput{s}{w}.*(abs(CInput{s}{w})>thresh*E{s}{w});  %大于阈值的保留
    end
end
denoise_S1=real(ifdct_wrapping(Ct,is_real,Nt,Nx));
%% %迭代去噪，输入为原始信号以及初始去噪结果
for i=1:5
        P=(Input-denoise_S1);
        inter=abs(P-median(P(:)));
        delta=median(inter(:))/0.675*1.345;
        E_out=siga(P,delta);
        Z=denoise_S1+E_out;   %得到该循环中奇异噪音衰减后的结果
               
        CInput=fdct_wrapping(Z,is_real,finest);     %进行曲波变换用的是实曲波,最外层为小波变换

        Smax=length(CInput);
        Sigma=sigma(i+1);     %求取噪音标准差，选取最大尺度

        Ct=CInput;
        for s=2:length(CInput)
            thresh=Sigma+Sigma*s;    
            for w=1:length(CInput{s})
                Ct{s}{w}=CInput{s}{w}.*(abs(CInput{s}{w})>thresh*E{s}{w});
            end
        end
        denoise_S1=real(ifdct_wrapping(Ct,is_real,Nt,Nx));
        i
end
%% 奇异值去除之后，残余随机噪音，最后再进行一次常规去噪
CInput=fdct_wrapping(denoise_S1,is_real,finest);    
Sigma=alpha*median(median(abs(CInput{Smax}{1})))/0.58*5;  

Ct=CInput;
for s=2:length(CInput)
    thresh=Sigma+Sigma*s;    
    for w=1:length(CInput{s})
        Ct{s}{w}=CInput{s}{w}.*(abs(CInput{s}{w})>thresh*E{s}{w}); 
    end
end
denoise_S1=real(ifdct_wrapping(Ct,is_real,Nt,Nx));
%% 输出去噪结果以及噪音信号
denoise_N=(Input-denoise_S1);
glwritesu(denoise_S1,0.001,'denoise_S.su');
glwritesu(denoise_N,0.001,'denoise_N.su');

figure;imagesc([Input,denoise_S1,denoise_N]*10^8);caxis([-1,1]);

