function [Ks, Cz] = C_Ks_parabol(Npart,L,dz,Ws,mKs)
%C_KS_LINEAR Analytical concentration proflie for a parabolitical turbidity profile
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
ustar = mKs/kappa*length(z)/sum(z.*(1-(z./L))); % total friction velocity (m.s⁻¹)
disp(['u* = ' num2str(ustar) 'm.s⁻¹'])

b = Ws/(kappa*ustar); % Rouse number 

za = 1; % ref depth (m)
Ca = Npart/(sum((z(2:end)/(L-z(2:end))).^(-b))*((L-za)/za)^(-b)*dz); % concentration at depth za (mps.m⁻³)

Ks = kappa*ustar*z.*(1-(z./L)); % linear Ks (m².s⁻¹)
Cz = Ca.*(z/za*(L-za)./(L-z)).^(-b); % concentration (mps.m⁻³)

subplot(1,2,1), plot(Ks,-z), xlabel('Ks (m².s⁻¹)'), ylabel('depth (m)')
subplot(1,2,2), plot(Cz,-z), xlabel('Cz (mps.m⁻³)'), ylabel('depth (m)')

end