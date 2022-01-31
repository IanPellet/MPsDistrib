function [K,dK] = Diffusivity(z,z_,dz,tension,cutoff,KZ_data,z0_data)
%DIFFUSIVITY Returns diffusivity profile
%   Kz data is interpolated on z_ and smoothed to avoid rappid changes in K
%   and dK. It returns the diffusivity K and its derivative allong the
%   depth z  at each points of z_.
%
% ARGUMENTS :
% z : double 1D array, water column disscretization, boundaries of the meshes (m)
% z_ : double 1D array, middle of the meshes (m)
% tension : double, tension parameter for the tension spline
% cutoff : double, cutt-off parameter for the tension spline
% KZ_data : double 1D array, K data to intepolate (m².s⁻¹)
% z0_data : double 1D array, depth of data points in KZ_data (m)
%
% RETURNS :
% K : double 1D array, diffusivity at each point of z_ (m².s⁻¹)
% dK : double 1D array, diffusivity gradiant at each point of z_ (m.s⁻¹)
%

% % If no data was passed : use data in DensiteFevrierRhoma
% if nargin == 5
%     DensiteFevrierRhoma
%     KZ_data = KZ_Fev10;
%     z0_data = z0;
% end

if -z(end) == z0_data(1)
    KZ_data(1) = 0;
else
    z0_data = reshape(z0_data, length(z0_data),1);
    z0_data = [-z(end) ; z0_data];
    KZ_data = [0 KZ_data];
end

Kinter = interp1(z0_data, KZ_data, -z,'linear');
Kinter(isnan(Kinter))=0;

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

% Compute diffusivity gradiant
dK = diff(Kspl)/dz;
% 
% figure(1), clf
% plot(K,-z_)
% figure(2), clf
% plot(dK,-z_)
end

