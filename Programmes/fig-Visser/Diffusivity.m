function [K,dK] = Diffusivity(z,z_,dz)
%DIFFUSIVITY Returns diffusivity profile
%   Pt Somlit, 10 Fev

addpath('../');

DensiteFevrierRhoma 

Kinter = interp1(z0, KZ_Fev10, -z,'linear');

%% 8-point moving average
npts = 8;
wts = ones(1,npts)*1/npts; % convolution matrix
Kmavg = conv(Kinter,wts,'valid');
zmavg = conv(z,wts,'valid');

Kspl = spline1d (z, zmavg, Kmavg, [min(z) max(z)], [0 0], 0.5)';
K = spline1d (z_, zmavg, Kmavg, [min(z_) max(z_)], [0 0], 0.5)';

dK = diff(Kspl)/dz;


end

