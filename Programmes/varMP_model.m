function [z_, CfinalNorm] = varMP_model(D, rhoP, nPart, dt, tf, dt_test, WindSpeed, month, Lon0, Lat0, path)
%VARMP_MODEL Summary of this function goes here
%   Detailed explanation goes here

fprintf(['\n\n\n--------------------- D = ' num2str(D)...
    ' -- rhoP = ' num2str(rhoP) ' ---------------------\n\n'])

clear Concentration err

 ModeleHydro='2012RHOMA_arome_003.nc';
    SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
    load(SauvegardeModeleHydro)

%% Speed computation formulas
Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];
indNom=3;

%% Water column parameters
[I0,J0]=ReperePoint(Lon,Lat,Lon0,Lat0);
L = H0(I0,J0);
% L = 50;
N = fix(L);  dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

%% Speed initialisation
[KZ_day,Row_day,z_day,z__day] = KsSalTemp(WindSpeed, month);

[K,dK] = wcp_interpolation(z_day,KZ_day,-z_); % Diffusivity
% K_val = Ks;
% K = ones(size(z))*K_val;
% dK = zeros(size(z));

g = 9.81 ; %m.s-1 (gravitational acceleration)
nuw = 1.1*10^-6; %m2.s-1 (kinematic viscosity of sea water)
rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 

S=rhoP./rhow;
D_=((g*(abs(S-1))/nuw^2).^(1/3))*D;
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
u=Ws; u(rhoP<rhow)=-Ws(rhoP<rhow);

%% Particules matrix definition
% part = [z1 z2 ... zn ; 
%         u1 u2 ... un ;
%         K1 K2 ... Kn ;
%        dK1 dK2...dK3]
z_part = linspace(0, L, nPart);
index = max(1, cast(z_part/dz, 'uint32'));
u_part = u(index);
K_part = K(index);
dK_part = dK(index);
part = [z_part ; u_part ; K_part ; dK_part];


%% Plot initial conditions
h_init = histogram(part(1,:), "BinEdges", z, 'Visible', 'off').Values;
CiNorm = h_init/dz*N/nPart;

%% Simulation
%nStep = tf/dt+1;
nStep = tf/dt_test+1;
if cast(nStep, 'uint32') ~= nStep
    nStep = cast(nStep+1, 'uint32');
end
T = ones(nStep,1);
ChNorm = ones(nStep,N);

z_past = part(1,:);
step = 1;
ChNorm(step,:) = [CiNorm];

t=0; OnContinue=true;
T(step,1) = t;
while OnContinue
    
    % Time update
    t=t+dt;
    
    % Particules update
    part(1,:) = Step_Lagrangien(part(1,:), part(2,:), part(3,:), part(4,:), dt);
    index = max(1, cast(part(1,:)/dz, 'uint32'));
    part(2,:) = u(index);
    part(3,:) = K(index);
    part(4,:) = dK(index);
    
    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 ) 
    %if true
        step = step+1;
        T(step,1) = t;
        z_present = part(1,:);
        
        % Computation of the concentration of MPs in each mesh
        h_past = histogram(z_past, "BinEdges", z, 'Visible', 'off').Values;
        h_present = histogram(z_present, "BinEdges", z, 'Visible', 'off').Values;
        CpastNorm = h_past/dz*N/nPart;
        CpresentNorm = h_present/dz*N/nPart;
        ChNorm(step,:) = CpresentNorm;     
                
        dC = max(abs(CpresentNorm - CpastNorm)/dt);
        
        if (t>=tf || dC/sum(CpresentNorm)*100 <= 0.005)
           OnContinue = false;
        end
        
        if mod(t,60)<dt/2
            disp([' Temps : ' num2str(t/tf*100) '% - Ecart : ' num2str(dC/sum(CpresentNorm)*100), '%'])
        end
        
        z_past = z_present;
    end
end
CfinalNorm = CpresentNorm;
%% Plot end profile
prof = figure(1);
plot(CfinalNorm,-z_, 'DisplayName',...
    join(['D = ',num2str(D*1e6),'µm ; rhop = ',num2str(rhoP),'kg.m⁻³']))
legend('Location','best')
ylim([-L 0])
xlim([0 10])
grid on
inv = z(end:-1:1);
yticks(-inv)
xlabel("Normed concentration (mps.m⁻³)")
ylabel("Depth (m)")
ttl = join(['D = ',num2str(D*1e6),'µm ; rhop = ',num2str(rhoP),'kg.m⁻³']);
title(ttl)
hold off

% path = '../../Ian/Results/varMP/';
prof_name = ['profile_D',num2str(D), '_rhop',num2str(rhoP),...
    '_nPart',num2str(nPart), '_dt', num2str(dt),...
    '_tf', num2str(tf), '_dtest', num2str(dt_test),'_wind', num2str(WindSpeed),...
    '_month', num2str(month)];
exportgraphics(prof,[path,prof_name,'.eps'],'ContentType','vector');
savefig(prof,[path,'fig/',prof_name,'.fig']);

end

