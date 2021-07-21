% Lagrangian transportation model 
% 
% 

% DensiteFevrierRhoma % load diffusivity data
% clearvars -except KZ_Fev10 z_Fev10,

date = datetime(2020,02,10);

%% Speed computation formulas
Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];
indNom = 3; % Ahrens's speed computation formula is chosen

%% Equilibrium test parameters
tf= 60*60*24*2; % maximum simulation time
dt_test = 60*60; % time interval between images

%% Load hydrodinamic model data
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro, 'H0', 'Lon', 'Lat')

% %% Water column parameters
% % Find depth of the column
% Lon0 = 5.29; Lat0 = 43.24; %point Somlit
% [I0,J0] = ReperePoint(Lon,Lat,Lon0,Lat0); % data indices corresponding to the location
% L = H0(I0,J0); % depth
% clear H0 Lon Lat,
L = 50;

% Column discretisation
N = 50; % number of meshes
dz = L/N; % depth of meshes
z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh

%% Particules initialisation
D=350e-6; % Particle diametre
rop = 950; % Particle density
nPart = 1e3; % number of particles initialisated in the column
zPart = linspace(0,L,nPart); % initial particle's positions (uniform repartition)

%% Use diffusivity data 10fev

[KZ_day,Row_day,z_day,z__day,Sal_day,Temp_day] = KsSalTemp2020(NaN, date);
% KZ_day = KZ_Fev10;
% z_day = z_Fev10;
% Interpolation and smoothing of the diffusivity data
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day);

%% Speed initialisation
g = 9.81 ; %m.s-1 (gravitational acceleration)
nuw = 1.1*10^-6; %m2.s-1 (kinematic viscosity of sea water)
% rhoW = getCTDrhow('10fev',z); % measured water density for 10fev
rhoW = interp1(-z__day,Row_day,z,'pchip'); % density of sea water

% compute particle's speed
S = rop./rhoW;     
D_ = ((g*(abs(S-1))/nuw^2).^(1/3))*D;
Ws = eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']); % absolute speed
u = Ws; 
u(rop<rhoW) = -Ws(rop<rhoW); % speed direction


%% Particules matrix definition
% part = [z1 z2 ... zn ; 
%         u1 u2 ... un ;
%         K1 K2 ... Kn ;
%        dK1 dK2...dK3]
index = max(1, fix(zPart/dz)); % find the mesh of each particle
u_part = u(index); % particle speed
K_part = K(index); % diffusivity
dK_part = dK(index); 
part = [zPart ; u_part ; K_part ; dK_part]; % particles data structure

%% Time step init
dt = 10;
ddK = diff(dK)/dz;
dt = min(dt, abs(min(1./ddK)/10)); % condition dt<<min(1/ddK) 
dt = min(dt, dz/max(abs(u))); % condition dt < dz/max|u|

%% Simulation
t=0; OnContinue=true;
figure(1)
while OnContinue
    
    % Time update
    t=t+dt;
    
    % Particules update
    part(1,:) = Step_Lagrangien(part(1,:), part(2,:), part(3,:), part(4,:), dt, L); % update position
    index = max(1, fix(part(1,:)/dz)); % index of the mesh of each particle
    % Update speed/diffusivity for each particle
    part(2,:) = u(index);
    part(3,:) = K(index);
    part(4,:) = dK(index);

    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 ) 
        if (t+dt>tf)
           OnContinue = false;
        end

        disp([' Temps : ' num2str(t/3600/24) 'j'])

        % Plot particle position
        plot(-part(1,:),'.')
        ylim([-L 0])
%         grid on
%         inv = z(end:-1:1);
%         yticks(-inv)
        
        pause(0) % pause code to visualise plot

    end
end

%% Results

pp = part(1,:); % final particle position
ppSort = sort(pp); % sort positions

topBound = dz; % top bound of the mesh
j = 1; % index of the mesh
count = zeros(N,1); % number of particle per mesh
% Count number of particle per mesh
for i=1:length(ppSort)
    if ppSort(i) > topBound
        topBound = topBound+dz;
        j = j+1;
    end
    count(j) = count(j)+1;
end, clear i,

Conc = count/dz; % compute concentration in each mesh

Ccalc = C_analytical(mean(Ws), mean(K), z_, nPart, L); % compute analytical solution

% plot results
figure(2), clf,
hold on
plot(Ccalc, -z_, 'DisplayName', 'Analytical solution')
plot(Conc, -z_, 'DisplayName', 'Model output')
hold off
legend('Location', 'best')
xlabel("Concentration (mps.m⁻³)")
ylabel("Depth (m)")


