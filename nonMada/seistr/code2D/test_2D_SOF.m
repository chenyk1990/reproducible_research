%% Demo for structure-oriented filtering algorithm of 2D data
%
%  Copyright (C) 2019 
%  Hang Wang, Yangkang Chen, and co-authors, 2019
%  
%  Reference
%  H. Wang, Y. Chen, O. Saad, W. Chen, Y. Oboue, L. Yang, S. Fomel, and Y. Chen, 2021, A Matlab code package for 2D/3D local slope estimation and structural filtering: in press.

%% Generate synthetic data
clear;clc;close all;

w=str_ricker(30,0.001,0.1);
t=zeros(300,1000);
sigma=300;A=100;B=200;
for i=1:size(t,2)
k=floor(-A*exp(-(i-size(t,2)/2).^2/sigma.^2)+B);
if k>0&&k<=size(t,1)
    t(k,i)=1;
end
end
for i=1:size(t,2)
data(:,i)=conv(t(:,i),w);
end
data=data(:,1:10:end);
data=data./max(max(data));
scnoi=(rand(size(data))*2-1)*0.2;
dn=data+scnoi;


% plot figures
lim1=-1;lim2=1;
x1=100;y1=100;dx=400;dy=500;

figure;imagesc(data);ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

figure;imagesc(dn);ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

%% Slope estimation
dtemp=dn*0;%dtemp is the preprocessed data
for i=1:size(dn,1)
    dtemp(i,:)=smooth(dn(i,:),5);
end

[dip]=str_dip2d(dtemp);

%plot figures
lim1=-1.2;lim2=2.7;
x1=100;y1=100;dx=400;dy=500;
figure;imagesc(dip);colormap(jet);ax = gca;
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);


%% Structural smoothing
r=2;
eps=0.01;
order=2;
% dn is the input noisy data, d1 is the output smoothed data
d1=str_pwsmooth_lop2d(dn,dip,r,order,eps);



%plot figures
lim1=-1;lim2=1;
x1=100;y1=100;dx=400;dy=500;

figure;imagesc(d1);ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

figure;imagesc(dn-d1);ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

figure;imagesc(dtemp);ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

figure;imagesc(dn-dtemp);ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

snrn=str_snr(data,dn);
snr1=str_snr(data,d1);
snrc=str_snr(data,dtemp);


