clc;clear;close all;

is_real=1;           % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
finest=2;            % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
alpha=1.0;           % 噪音标准差的alpha倍阈值（1.2左右较为理想）

niter=10;

infile='~/chenyk.data2/various/cyksmall/zq_real_input.su';
[ori]=glreadsu(infile,'su');
dn=scalecyk(ori,2);
[n1,n2]=size(dn);
figure;imagesc([dn]);colormap(seis);caxis([-0.5,0.5]);


F=ones(n1,n2);                                  % ones(n)返回n*n的1矩阵，频率域
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

Cdn=fdct_wrapping(dn,is_real,finest);     %进行曲波变换用的是实曲波,最外层为小波变换
Smax=length(Cdn);
Sigma0=alpha*median(median(abs(Cdn{Smax}{1})))/0.58;     %求取噪音标准差，选取最大尺度
Sigma=Sigma0;
% sigma=[Sigma,5*Sigma,2*Sigma, Sigma, 0.6*Sigma,Sigma/5];
% sigma=[Sigma,5*Sigma,2*Sigma, Sigma, 0.6*Sigma,Sigma*0.1];
sigma=[Sigma,linspace(2.5*Sigma,0.5*Sigma,niter)];
Sigma=sigma(1);

Ct=Cdn;
for s=2:length(Cdn)
    thresh=Sigma+Sigma*s;    %最外层设置为4*sigma
    for w=1:length(Cdn{s})
        Ct{s}{w}=Cdn{s}{w}.*(abs(Cdn{s}{w})>thresh*E{s}{w});  %大于阈值的保留
    end
end
d1=real(ifdct_wrapping(Ct,is_real,n1,n2));

d2=d1;

% figure;imagesc(dn);colormap(seis);caxis([-0.5,0.5]);colormap(seis);
figure;imagesc([dn,d1,dn-d1]);colormap(seis);caxis([-0.5,0.5]);colormap(seis);

%% %迭代去噪，输入为原始信号以及初始去噪结果
Sigma=Sigma0;
% sigma=[Sigma,5*Sigma,2*Sigma, Sigma, 0.6*Sigma,Sigma/5];
% sigma=[Sigma,5*Sigma,2*Sigma, Sigma, 0.6*Sigma,Sigma*0.1];
sigma=[Sigma,linspace(2.5*Sigma,0.5*Sigma,niter)];


for i=1:niter
    P=(dn-d1);
    inter=abs(P-median(P(:)));
    delta=median(inter(:))/0.675*1.345
    E_out=siga(P,delta);
    Z=d1+E_out;   %得到该循环中奇异噪音衰减后的结果
    
    Cdn=fdct_wrapping(Z,is_real,finest);     %进行曲波变换用的是实曲波,最外层为小波变换
    
    Smax=length(Cdn);
    Sigma=sigma(i+1);     %求取噪音标准差，选取最大尺度
    
    Ct=Cdn;
    for s=2:length(Cdn)
        thresh=Sigma+Sigma*s;
        for w=1:length(Cdn{s})
            Ct{s}{w}=Cdn{s}{w}.*(abs(Cdn{s}{w})>thresh*E{s}{w});
        end
    end
    d1=real(ifdct_wrapping(Ct,is_real,n1,n2));
    
    pause(1);figure(2);imagesc([dn,d1,dn-d1]);colormap(seis);caxis([-0.5,0.5]);
end



%% SOSVMF
dipn=dip2d_shan(d1,5,20,2,0.01,[20,5,1]);
figure;imagesc([dipn,dipn]);colorbar;colormap(jet);caxis([-1,2]);

dt=0.004;
t=[0:n1-1]*dt;x=[1:n2];

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dipn);colorbar;colormap(jet);caxis([-1,2]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_dip_dn.eps');



type_mf=1;ifsmooth=0;
ns=3;
[~,d3]=pwsmooth_lop_mf(0,0,dipn,[],n1,n2,ns,2,0.01,n1*n2,n1*n2,type_mf,ifsmooth,dn,[]);
d3=reshape(d3,n1,n2);


figure;imagesc([dn,d1,dn-d1]);colormap(seis);caxis([-0.5,0.5]);colormap(seis);%new
figure;imagesc([dn,d2,dn-d2]);colormap(seis);caxis([-0.5,0.5]);colormap(seis);%curvelet
figure;imagesc([dn,d3,dn-d3]);colormap(seis);caxis([-0.5,0.5]);colormap(seis);%curvelet


%% %迭代去噪，输入为原始信号以及初始去噪结果
d4=d3;
nfw=[3,3,3,3,3];
niter=10;
nfw=3*ones(niter,1);
% Sigma=0.2104*0.5;
Sigma=Sigma0;
% sigma=[Sigma,linspace(2.5*Sigma,0.25*Sigma,niter)];%max 10.02
sigma=[Sigma,linspace(Sigma,0.1*Sigma,niter)];%very good performance
% snrs4=zeros(niter+1,1);
snrs4=[];

% dipn=dip2d_shan(dn,5,20,2,0.01,[20,5,1]);

for i=1:niter-2
    P=(dn-d4);
    inter=abs(P-median(P(:)));
    delta=median(inter(:))/0.675*1.345
    E_out=siga(P,delta);
    Z=d4+E_out;   %得到该循环中奇异噪音衰减后的结果
    
    [~,Z]=pwsmooth_lop_mf(0,0,dipn,[],n1,n2,nfw(i),2,0.01,n1*n2,n1*n2,type_mf,ifsmooth,dn,[]);
    Z=reshape(Z,n1,n2);
    
    Cdn=fdct_wrapping(Z,is_real,finest);     %进行曲波变换用的是实曲波,最外层为小波变换
    
    Smax=length(Cdn);
    Sigma=sigma(i+1);     %求取噪音标准差，选取最大尺度
    
    Ct=Cdn;
    for s=2:length(Cdn)
        thresh=Sigma+Sigma*s;
        for w=1:length(Cdn{s})
            Ct{s}{w}=Cdn{s}{w}.*(abs(Cdn{s}{w})>thresh*E{s}{w});
        end
    end
    
    
%     dipn=dip2d_shan(d4,5,20,2,0.01,[20,5,1]);
    d4=real(ifdct_wrapping(Ct,is_real,n1,n2));


    pause(1);figure(2);imagesc([dn,d4,dn-d4]);colormap(seis);caxis([-0.5,0.5]);

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dipn);colorbar;colormap(jet);caxis([-1,2]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
tname=strcat(['f_dip_',num2str(i),'.eps']);
print(gcf,'-depsc','-r200',tname);

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,d4);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
tname=strcat(['f_d_',num2str(i),'.eps']);
print(gcf,'-depsc','-r200',tname);

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dn-d4);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
tname=strcat(['f_n_',num2str(i),'.eps']);
print(gcf,'-depsc','-r200',tname);
end
figure;plot(snrs4);



%% plot

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dn);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_dn.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,d2);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_curv.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,d1);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_curvi.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,d3);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_sosvmf.eps');

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,d4);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_sosvmfi.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dn-d2);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_curv_n.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dn-d1);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_curvi_n.eps');


figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dn-d3);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_sosvmf_n.eps');

figure('units','normalized','Position',[0.2 0.4 0.4, 0.8],'color','w');
imagesc(x,t,dn-d4);colormap(seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_sosvmfi_n.eps');



rect=[10,5,1];niter=20;eps=0;verb=0;
[simi11]=localsimi(dn-d1,d1,rect,niter,eps,verb);
[simi22]=localsimi(dn-d2,d2,rect,niter,eps,verb);
[simi33]=localsimi(dn-d3,d3,rect,niter,eps,verb);
[simi44]=localsimi(dn-d4,d4,rect,niter,eps,verb);

figure('units','normalized','Position',[0.0 0.0 0.6 1.0],'color','w');
imagesc(x,t,simi11);colormap(jet);
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 30;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_simi2_curvi.eps');


figure('units','normalized','Position',[0.0 0.0 0.6 1.0],'color','w');
imagesc(x,t,simi22);colormap(jet);
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 30;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_simi2_curv.eps');


figure('units','normalized','Position',[0.0 0.0 0.6 1.0],'color','w');
imagesc(x,t,simi33);colormap(jet);
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 30;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_simi2_sosvmf.eps');


figure('units','normalized','Position',[0.0 0.0 0.6 1.0],'color','w');
imagesc(x,t,simi44);colormap(jet);
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 30;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',30);
xlabel('Trace','Fontsize',30);
% title('Noise','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
print(gcf,'-depsc','-r200','f_simi2_sosvmfi.eps');




