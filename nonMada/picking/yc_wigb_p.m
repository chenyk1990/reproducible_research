function yc_wigb_p(a,pk,style,scal,x,z,amx)
%YC_WIGB_P: Plot seismic data using wiggles with arrival picker
%
%  WIGB(a,pk,scal,x,z,amx,style)
%
%  IN  pk:arrival picker
%       a:     seismic data (a matrix, traces are columns)
%      scale: multiple data by scal
%      x:     horizontal axis (often offset)
%      z:     vertical axis (often time)
%
%  Note
%
%    If only 'a' is enter, 'scal,x,z,amn,amx' are set automatically; 
%    otherwise, 'scal' is a scalar; 'x, z' are vectors for annotation in 
%    offset and time, amx are the amplitude range.
%
%
%  Author(s): Xingong Li (Intseis, Integrated Seismic Solutions)
%  Copyright 1998-2003 Xingong 
% 
%  Modified by Yangkang Chen, Nov, 2016
% 
% 
% Reference:
% Chen, Y., 2020, Automatic microseismic event picking via unsupervised
% machine learning, GJI, 1750–1764
%
% Related reference:
% Chen, Y., 2018, Fast waveform detection for microseismic imaging using unsupervised machine learning, GJI, 1185-1199.
% Zhang et al., 2020, Convolutional Neural Networks for Microseismic Waveform Classification and Arrival Picking, Geophysics, WA227–WA240.
% Chen, et al., 2019, Automatic Waveform Classification and Arrival Picking Based on Convolutional Neural Network, ESS, 1244-1261.
% Saad and Chen, 2020, Automatic waveform-based source-location imaging using deep learning extracted microseismic signals, Geophysics, KS171–KS183.
% Qu et al., 2020, Automatic high-resolution microseismic event detection via supervised machine learning, GJI, 1881–1895.
%
% For earthquake detection and picking
% Saad and Chen, 2021, CapsPhase: Capsule Neural Network for Seismic Phase Classification and Picking, TGRS.
% Saad et al., 2021, SCALODEEP: A Highly Generalized Deep Learning Framework for Real-time Earthquake Detection, JGR, e2020JB021473
% Saad and Chen, 2020, Earthquake Detection and P-wave Arrival Time Picking using Capsule Neural Network, TGRS.


 if nargin == 0, nx=10;nz=10; a = rand(nz,nx)-0.5; end;

 [nz,nx]=size(a);

 trmx= max(abs(a));
 if (nargin <= 6);   amx=mean(trmx);end;
 if (nargin <= 4);   x=[1:nx]; z=[1:nz];  end;
 if (nargin <= 3); scal =1;end;
 if (nargin <= 2); style='ro'; end;

 if (nargin <= 1); error('Missing parameters!'); end

 
 if nx <= 1; disp(' ERR:PlotWig: nx has to be more than 1');return;end;

 % take the average as dx

	dx1 = abs(x(2:nx)-x(1:nx-1));
 	dx = median(dx1);


 dz=z(2)-z(1);
 xmx=max(max(a)); xmn=min(min(a)); 

 if scal == 0; scal=1; end;
 a = a * dx /amx; 
 a = a * scal;

 fprintf(' PlotWig: data range [%f, %f], plotted max %f \n',xmn,xmx,amx);
 
% set display range 

 x1=min(x)-2.0*dx; x2=max(x)+2.0*dx;
 z1=min(z)-dz; z2=max(z)+dz;
 
 set(gca,'NextPlot','add','Box','on', ...
  'XLim', [x1 x2], ...
  'YDir','reverse', ...
  'YLim',[z1 z2]);

 

	fillcolor = [0 0 0];
	linecolor = [0 0 0];
	linewidth = 1.;

	z=z'; 	% input as row vector
	zstart=z(1);
	zend  =z(nz);

for i=1:nx,
   
  if trmx(i) ~= 0;    % skip the zero traces
	tr=a(:,i); 	% --- one scale for all section
  	s = sign(tr) ;
  	i1= find( s(1:nz-1) ~= s(2:nz) );	% zero crossing points
	npos = length(i1);


	%12/7/97 
	zadd = i1 + tr(i1) ./ (tr(i1) - tr(i1+1)); %locations with 0 amplitudes
	aadd = zeros(size(zadd));

	[zpos,vpos] = find(tr >0);
	[zz,iz] = sort([zpos; zadd]); 	% indices of zero point plus positives
	aa = [tr(zpos); aadd];
	aa = aa(iz);

	% be careful at the ends
		if tr(1)>0, 	a0=0; z0=1.00;
		else, 		a0=0; z0=zadd(1);
		end;
		if tr(nz)>0, 	a1=0; z1=nz; 
		else, 		a1=0; z1=max(zadd);
		end;
			
	zz = [z0; zz; z1; z0];
 	aa = [a0; aa; a1; a0];
		
	zzz = zstart + zz*dz -dz;

patch( aa+x(i) , zzz,  fillcolor);

	line( 'Color',[1 1 1],  ...
         'Xdata', x(i)+[0 0], 'Ydata',[zstart zend]); % remove zero line

%'LineWidth',linewidth, ...
%12/7/97 'Xdata', x(i)+[0 0], 'Ydata',[z0 z1]*dz);	% remove zero line

	line( 'Color',linecolor,  ...
	 'LineWidth',linewidth, ...
	 'Xdata', tr+x(i), 'Ydata',z);	% negatives line
 z0=z(1);
 dz=z(2)-z(1);
    plot(x(i),z0+dz*(pk(i)-1),style,'linewidth',2);
    
   else % zeros trace
line( 'Color',linecolor,  ...
	 'LineWidth',linewidth, ...
         'Xdata', [x(i) x(i)], 'Ydata',[zstart zend]);
   end;
end;
