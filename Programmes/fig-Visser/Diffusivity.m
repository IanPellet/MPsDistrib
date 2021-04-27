function [K_,dK] = Diffusivity(z,z_,dz)
%DIFFUSIVITY Returns diffusivity profile
%   Pt Somlit, 10 Fev
% H : depth of the water column (m)
% N : number of meshes

addpath('../');

DensiteFevrierRhoma 
% [K,dK] = wcp_interpolation(z0,KZ_Fev10,-z); % Diffusivity
K_ = interp1(z0, KZ_Fev10, -z_,'linear');
K = interp1(z0, KZ_Fev10, -z,'linear');

%disp(size(Kz))

% %% 8-point moving average
% npts = 8;
% wts = ones(1,npts)*1/npts; % convolution matrix
% Kz = [zeros(1,npts/2) Kz zeros(1,npts/2+1)]; % padding
% Kz = conv(Kz,wts,'valid');
% %disp(size(Kz))

dK = diff(K)/dz;
%disp(size(dKz))
%dKz = conv(dKz,wts,'valid');

end

