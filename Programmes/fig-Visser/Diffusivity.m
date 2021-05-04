function [K,dK] = Diffusivity(z,z_,dz,tension,cutoff,KZ_Fev10,z0)
%DIFFUSIVITY Returns diffusivity profile
%   Pt Somlit, 10 Fev

if nargin == 5
    addpath('../');
    DensiteFevrierRhoma
end


z0 = [-z(end) ; z0];
KZ_Fev10 = [0 KZ_Fev10];

Kinter = interp1(z0, KZ_Fev10, -z,'linear');

%% 8-point moving average
npts = 8;
wts = ones(1,npts)*1/npts; % convolution matrix
Kmavg = conv(Kinter,wts,'valid');
zmavg = conv(z,wts,'valid');

% tension = 0.5;
% cutoff = 0.1;
Kspl = spline1d (z, zmavg, Kmavg, [min(z) max(z)], [0 0], tension, cutoff)';
K = spline1d (z_, zmavg, Kmavg, [min(z_) max(z_)], [0 0], tension, cutoff)';
Kspl(Kspl<0) = 0;
K(K<0) = 0;

dK = diff(Kspl)/dz;
% 
% figure(1), clf
% plot(K,-z_)
% figure(2), clf
% plot(dK,-z_)
end

