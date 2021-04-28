function [K_,dK] = Diffusivity(z,z_,dz)
%DIFFUSIVITY Returns diffusivity profile
%   Pt Somlit, 10 Fev
% H : depth of the water column (m)
% N : number of meshes

addpath('../');

DensiteFevrierRhoma 
% [K,dK] = wcp_interpolation(z0,KZ_Fev10,-z); % Diffusivity
% K_ = interp1(z0, KZ_Fev10, -z_,'spline');
K = interp1(z0, KZ_Fev10, -z,'spline');
%% 2-point moving average
npts = 2;
wts = ones(1,npts)*1/npts; % convolution matrix
% disp(wts)
K_ = conv(K,wts,'valid');



% K_(1:5) = K_(5);
% K(1:5) = K(5);

% disp(size(K_))

%% 9-point moving average
% npts = 51;
% wts = ones(1,npts)*1/npts; % convolution matrix
% disp(wts)
% K1_ = [zeros(1,(npts-1)/2) K_ zeros(1,(npts-1)/2)]; % padding
% K2_ = conv(K1_,wts,'valid');
% K1 = [zeros(1,(npts-1)/2) K zeros(1,(npts-1)/2)]; % padding
% K2 = conv(K1,wts,'valid');


dK = diff(K)/dz;
% dK2 = diff(K2)/dz;
% dKz = conv(dKz,wts,'valid');
% 
% figure(2), 
% ax1 = subplot(1,2,1);
% plot(ax1,K,-z,'b',K_,-z_,'m',KZ_Fev10, z0,'xg')
% ax2 = subplot(1,2,2);
% plot(ax2,dK,-z_)

end

