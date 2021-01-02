function  yc_wigbh(din,x,z,clip)
%YC_WIGBH wigb horizontally
%din: size nz*nx
%
% Nature of clip: the plotting amplitude correspoding to the clip value is 1 
% 
%Written by by Yangkang Chen, Nov, 2018
% 
% Example:
% d=levents;figure;wigbh(d,1:50,0.004*[0:500],2);


[nz,nx]=size(din);
vmax=max(din(:));
if nargin==1
   x=1:nx;
   z=1:nz; 
   clip=vmax;
end
if nargin==3
   clip=vmax;
end

din=din*vmax/clip;


for i = 1 : nx
    plot(z, din(:,i)+ x(i),'k','LineWidth',1.5); hold on
end
% set(gca,'ytick',(25:5:95)*2);
% set(gca,'yticklabel',[25:5:95]);
set(gca,'FontSize',20);
%set(gca,'YDir','reverse');
set(gca,'linewidth',1.5);
axis tight
% xlabel('time/s');
% xlim([20 60]);
% ylabel('gcarc');
set(gcf,'color','white');

end

