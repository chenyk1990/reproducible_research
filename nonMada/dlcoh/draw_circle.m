function h = draw_circle(x,y,r,c)%c: color
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,c,'linewidth',4);
hold off
