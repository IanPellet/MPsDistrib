function [Kz, dKz] = wcp_interpolation(z_mes, Kz_mes, z)
%WCP_INTERPOLATION Water Column Profile interpolation
%   Data points are first proliferated with [linear] interpolation then
%   smoothed with [8-point] moving average

%% Linear interpolation 
Kz = interp1(z_mes, Kz_mes, z,'linear');

disp(size(Kz))

%% 8-point moving average
npts = 8;
wts = ones(1,npts); % convolution matrix
Kz = [zeros(1,npts/2) Kz zeros(1,npts/2+1)]; % padding
Kz = conv(Kz,wts,'valid');
disp(size(Kz))

dKz = diff(Kz);
disp(size(dKz))
%dKz = conv(dKz,wts,'valid');

end

