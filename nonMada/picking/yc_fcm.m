function [ center, U, objs ] = yc_fcm( data, Nc, params)
% 
% INPUT:
% data:  input data (N*nf) (feacture vectors: nf features)   
% Nc:       number of clusters
% params:   parameter vector (same as FCM)
%   
%       params(1): exponent for the matrix U             (default: 2.0)
%       params(2): maximum number of iterations          (default: 100)
%       params(3): minimum amount of improvement         (default: 1e-5)
%       params(4): info display during iteration         (default: 1)
%
% OUTPUT:
% center:   coordinates for each cluster center (each row for each cluster)
% U: The membership function matrix U contains the grade of membership of
%   each DATA point in each cluster. Grades between 0 and 1 indicate that the
%   data point has partial membership in a cluster (potential: 1-> high potential the point falls into the cluster).
% obj: objective function value during iterations
% 
% Copyright: Yangkang Chen, July, 2017
% 
% Example: see test/test_micro_fcm_cyk.m
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

if nargin == 2
params= [2;	% exponent for the partition matrix U (m)
		100;	% max. number of iteration
		1e-5;	% min. amount of improvement
		1];	% info display during iteration 
end


[N,Nf]=size(data);
%N: number of data points
%Nf: number of features (feature vectors)

m = params(1);		% Exponent for U 
niter= params(2);		% Max. iteration
eps = params(3);		% Min. improvement (convergence paramter)
verb = params(4);		% verbosity or not

objs= zeros(niter, 1);

%% initialize membership matrix U (randomize and normalize)
rand('state',201718);
U = rand(Nc, N);
col_sum = sum(U);
U = U./col_sum(ones(Nc, 1), :); % spread col_sum to the same size as U (sfspray axis=1 n=Nc)

%% 
for iter = 1:niter,

%% update membership matrix
%U,mf: Nc*N
%data: N*Nf
%center:Nc*Nf

mf = U.^m;       % MF matrix after exponential modification (mf: Nf*N)
center = mf*data./((ones(size(data, 2), 1)*sum(mf'))'); 
% new center (denominator is to make matrix same size as the nominator: Nf*Nf)

dist = distfcm(center, data);       % fill the distance matrix

dist=zeros(Nc,N);

%% calcualate distance
if Nf > 1,
    for k = 1:Nc
	dist(k, :) = sqrt(sum(((data-ones(N, 1)*center(k, :)).^2)'));
    end
else	% 1-D data (Nf=1, for 1D data clustering)
    for k = 1:Nc,
	dist(k, :) = abs(center(k)-data)';
    end
end

objs(iter) = sum(sum((dist.^2).*mf));  % objective function (L2 norm based distance)

tmp = dist.^(-2/(m-1));      % calculate new U, suppose m != 1
U = tmp./(ones(Nc, 1)*sum(tmp)); %sfspray axis=1 n=Nc

	if verb, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', iter, objs(iter));
	end
	
	% check termination condition
	if iter > 1,
		if abs(objs(iter) - objs(iter-1)) < eps, break; end,
	end
end

return























end

