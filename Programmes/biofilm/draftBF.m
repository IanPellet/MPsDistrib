D = 350e-6; % Particle diameter (m)

dt = 10;
tf = 60*60*24*10;
t = 0:dt:tf;
tday = t/60/60/24;

H = 70;
dz = 1;
z = 0:dz:H;


%% Vertical chlorophyll-a profile
avgChlaZbase = [0.091 0.151 0.185 0.250 0.338 0.410 0.578 1.206 2.950];
Zbase = [119.1 99.9 91 80.2 70.3 63.4 54.4 39.8 26.1]; % Euphotic zone depth
Cb = [0.471 0.533 0.428 0.570 0.611 0.390 0.569 0.835 0.188 ]; % Normalized surface value
s = [0.135 0.172 0.002 0.173 0.214 0.109 0.183 0.298 0.000]; % Normalized slope
Cmax = [1.572 1.194 1.015 0.766 0.676 0.788 0.608 0.382 0.885]; % Normalized peak concentration
Zmax = [0.969 0.921 0.905 0.814 0.663 0.521 0.452 0.512 0.378]; % Normalized depth of the peak
deltaZ = [0.393 0.435 0.630 0.586 0.539 0.681 0.744 0.625 1.081]; % Normalized width of the peak

zn = z./H;
% for i = 1:length(s)
%     Chlz = Cb(i) - s(i).*zn + Cmax(i).*exp(-((zn-Zmax(i))./deltaZ(i)).^2);
%     plot(Chlz, -z, 'DisplayName', num2str(i))
%     hold on
% end
% hold off
% legend('Location', 'best')
i = 6;
Chlz = Cb(i) - s(i).*zn + Cmax(i).*exp(-((zn-Zmax(i))./deltaZ(i)).^2);
Chlaz = Chlz * avgChlaZbase(i); % mg.m⁻³

% Chla_C = 0.003 + 1.0154 * exp(0.050*T) * exp(0.059*Iz/1e6) * muPrime; 

%% Light intensity profile
Im = 1.2e8 ; % Surface light intensity at noon (μE.m⁻².d⁻¹) 
% Surface light intensity (μE.m⁻².d⁻¹) 
I0 = Im*sinpi(2*tday); % 12h of day light
I0(I0<0) = 0; % 12h of night with no light

kw = 0.2;
kp = 0.02;

kz = kw + kp.*Chlaz; % light extinction

Iz = exp(-kz.*z)' .* ones(length(z),1)*I0;

%% Temperature profile
Tsurf = 25; % Temperature at ocean surface (°C)
Tbot = 1.5; % Temperature at ocean bottom (°C)
p = 2; % Steepness temperature decrease
zc = -300; % Depth thermocline (m)
Tz = Tsurf + (Tbot-Tsurf).*z.^p./(z.^p+zc.^p);

%% Algae growth
muMax = 1.85; % Maximum growth rate algae (d⁻¹)
Iopt = 1.75392e+13; % Optimal light intensity algae growth (μE.m⁻².d⁻¹)
alpha = 0.12; % Initial slope (d⁻¹)

muOptI = muMax .* Iz ./ (Iz+muMax/alpha.*(Iz./Iopt-1).^2);

Tmin = 0.2; % Minimum temperature algae growth (°C)
Topt = 26.7; % Optimal temperature algae growth (°C)
Tmax = 33.3; % Maximum temperature algae growth (°C)

PhiT = ((Tz-Tmax).*(Tz-Tmin).^2) ./ ((Topt-Tmin) .* ((Topt-Tmin) .* (Tz-Topt)...
    - (Topt-Tmax) .* (Topt+Tmin-2.*Tz)));

muTI = muOptI(:,6).*PhiT';
muTI(Tz > Tmax | Tz < Tmin) = 0;

%% Fouling through collision
betaAbrownian = 4*pi*(Dpl+DA) * (rtot+rA);
betaAsetling = 1.5*pi*rA^2*Vs;
betaAshear = 1.3*gamma*(rtot+rA)^3;
betaA = betaAbrownian + betaAsetling + betaAshear;

omegaPl = 4 * pi * (D/2)^2;