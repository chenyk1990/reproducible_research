%% Demo for structure-oriented filtering algorithm of 3D data
%
%  Copyright (C) 2019 
%  Hang Wang, Yangkang Chen, and co-authors, 2019
%  
%  Reference
%  H. Wang, Y. Chen, O. Saad, W. Chen, Y. Oboue, L. Yang, S. Fomel, and Y. Chen, 2021, A Matlab code package for 2D/3D local slope estimation and structural filtering: in press.

%% load 3D data
clear;clc;close;

load real3d.mat;

% plot 
lim1=-1;lim2=1;
x1=100;y1=100;dx=400;dy=500;
figure;imagesc(cmp(:,:,5));ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

x=permute(cmp,[3,2,1]);
x1=100;y1=100;dx=500;dy=500;
figure;h=slice(x,25,5,50,'cubic');%
ax = gca;load('MyColormaps','mycmap');
colormap(ax,mycmap);set(h,'FaceColor','interp','EdgeColor','none');
set(ax, 'CLim', [lim1 lim2]);set(gca,'ZDir','reverse');
set(gcf,'position',[x1,y1,dx,dy]);colorbar;
xlabel('Xline','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Inline','FontName','Arial','FontWeight','Bold','FontSize',14);
zlabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);


cmpn=cmp;
%% 3D slope calculation (inline and xline)
[dipi,dipx] = str_dip3d(cmpn);

%plot figures
x1=100;y1=100;dx=400;dy=500;
figure;imagesc(dipi(:,:,5));colormap(jet);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

%% Structural smoothing
r1=2;
r2=2;
eps=0.01;
order=2;

cmpn_d1=str_pwsmooth_lop3d(cmpn,dipi,dipx,r1,r2,eps,order);



%plot figures
lim1=-1;lim2=1;
x1=100;y1=100;dx=400;dy=500;

figure;imagesc(cmpn_d1(:,:,5));ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

figure;imagesc(cmpn(:,:,5)-cmpn_d1(:,:,5));ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

temp=cmpn*0;%temp is the smoothed data by the conventional method
for i=1:size(cmpn,1)
    for j=1:size(cmpn,2)
        temp(i,j,:)=smooth(cmpn(i,j,:));
    end
end

for i=1:size(cmpn,1)
    for j=1:size(cmpn,3)
        temp(i,:,j)=smooth(temp(i,:,j));
    end
end

figure;imagesc(temp(:,:,5));ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

figure;imagesc(cmpn(:,:,5)-temp(:,:,5));ax = gca;
load('MyColormaps','mycmap');colormap(ax,mycmap);
set(ax, 'CLim', [lim1 lim2]);
set(gcf,'position',[x1,y1,dx,dy]);
colorbar;xlabel('Trace','FontName','Arial','FontWeight','Bold','FontSize',14);
ylabel('Time (ms)','FontName','Arial','FontWeight','Bold','FontSize',14);
set(gca,'FontName','Arial','FontSize',14,'LineWidth',1);

