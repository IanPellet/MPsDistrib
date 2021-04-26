function [Ks, Cz] = C_Ks_linear(Npart,L,dz,Ws,mKs)
%C_KS_LINEAR Analytical concentration proflie for a linear turbidity profile
% Npart : number of particles
% L : depth of the water column (m)
% dz : meshes depth (m)
% Ws : settling velocity (m.s⁻¹)
% mKs :  mean Ks (m².s⁻¹)


% Npart = 50e3; % number of particles
% dz = 1; % meshes depth (m)
% Ws = 0.001; % settling velocity (m.s⁻¹)
% mKs = 1; % mean Ks (m².s⁻¹)

z = 0:dz:L; % depth (m)

kappa = 0.40; % Von Karman's constant
ustar = mKs/kappa*length(z)/sum(z); % total friction velocity (m.s⁻¹)
disp(['u* = ' num2str(ustar) 'm.s⁻¹'])

b = Ws/(kappa*ustar); % Rouse number 

za = 1; % ref depth (m)
Ca = Npart/(sum(z(2:end).^(-b))*dz)*za^(-b); % concentration at depth za (mps.m⁻³)

Ks = kappa*ustar*z; % linear Ks (m².s⁻¹)
Cz = Ca*(z/za).^(-b); % concentration (mps.m⁻³)

subplot(1,2,1), plot(Ks,-z), xlabel('Ks (m².s⁻¹)'), ylabel('depth (m)')
subplot(1,2,2), plot(Cz,-z), xlabel('Cz (mps.m⁻³)'), ylabel('depth (m)')

end