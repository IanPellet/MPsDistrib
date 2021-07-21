function [z_, CfinalNorm, PartPos, ppHistory] = varMP_model(D, rhoP, type, nPart, tf, dt_test, WindSpeed, month, Lon0, Lat0, L, N, day, saveNsec)
%VARMP_MODEL Summary of this function goes here
%   Detailed explanation goes here

global g nuw rhow 

fprintf(['\n\n\n--------------------- D = ' num2str(D*1e6) 'µm'...
    ' -- rhoP = ' num2str(rhoP) 'kg.m⁻³ -- Type = ' num2str(type) ' ---------------------\n\n'])

rng('shuffle'); % seed init based on clock

clear Concentration err

 ModeleHydro='2012RHOMA_arome_003.nc';
    SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
    load(SauvegardeModeleHydro)

%% Speed computation formulas
Nom=[... 
    ;{'Khatmullina'}...%Khatmullina (2016)
    ;{'Ahrens'}...%Ahrens (2000)
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
   ];
indNom = 2;
% indNom = type+1;
% if type > 1
%     disp("WARNING : Speed formula not implemented for this type of particle !");
% end

%% Water column parameters
% [I0,J0]=ReperePoint(Lon,Lat,Lon0,Lat0);
% L = H0(I0,J0);
% disp(L);
% L = 56;
% N = fix(L);  
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

%% Speed initialisation
% [KZ_day,Row_day,z_day,z__day] = KsSalTemp(WindSpeed, month);

DensiteFevrierRhoma
if nargin < 4
    day = '10fev';
end
if  strcmp(day, '10fev')
    KZ_day = KZ_Fev10;
    z_day = z_Fev10;
elseif  strcmp(day, '3fev')
    KZ_day = KZ_Fev03;
    z_day = z0;
else
    error('Day argument not recognised');
    %     [KZ_day,Row_day,z_day,z__day] = KsSalTemp(WindSpeed, month);
end

addpath('./fig-Visser/');
% [K,dK] = wcp_interpolation(z_day,KZ_day,-z_); % Diffusivity
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day);


g = 9.81 ; %m.s-1 (gravitational acceleration)
nuw = 1.1*10^-6; %m2.s-1 (kinematic viscosity of sea water)
rhow = getCTDrhow(day,z);

g_red = g.*abs(rhoP-rhow)./rhow;
S=rhoP./rhow;
D_=((g*(abs(S-1))/nuw^2).^(1/3))*D;
l = 0.5e-3;
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_,l,g_red);']);
u=Ws; u(rhoP<rhow)=-Ws(rhoP<rhow);

%% Time step init
dt = 10;
ddK = diff(dK)/dz;
dt = min(dt, abs(min(1./ddK)/10)); % condition dt<<min(1/ddK) 
dt = min(dt, dz/max(abs(u))); % condition dt < dz/max|u|

%% Particules matrix definition
% part = [z1 z2 ... zn ; 
%         u1 u2 ... un ;
%         K1 K2 ... Kn ;
%        dK1 dK2...dK3]
z_part = linspace(0, L, nPart);
index = max(1, fix(z_part/dz));
u_part = u(index);
K_part = K(index);
dK_part = dK(index);
part = [z_part ; u_part ; K_part ; dK_part];


%% Initial conditions
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

if nargin > 13 && saveNsec ~= 0
    saveNstep = fix(saveNsec/dt)+1;
    ppHistory = NaN(saveNstep,nPart);
    saveStep = 0;
end

t=0; OnContinue=true;
T(step,1) = t;
while OnContinue
    
    % Time update
    t=t+dt;
    
    % Particules update
    part(1,:) = Step_Lagrangien(part(1,:), part(2,:), part(3,:), part(4,:), dt, L);
    index = max(1, cast(part(1,:)/dz, 'uint32'));
    part(2,:) = u(index);
    part(3,:) = K(index);
    part(4,:) = dK(index);
    
    if t >= tf-saveNsec
        saveStep = saveStep+1;
        ppHistory(saveStep,:) = part(1,:);
    end
    
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
        
        if (t+dt>=tf || dC/sum(CpresentNorm)*100 <= 0.0005)
           OnContinue = false;
        end
        
        if mod(t,60)<dt/2
            disp([' Temps : ' num2str(t/tf*100) '% - Ecart : ' num2str(dC/sum(CpresentNorm)*100), '%'])
        end
        
        z_past = z_present;
    end
end
CfinalNorm = CpresentNorm;

PartPos = part(1,:);

%% Plot end profile
% figure(1)
% plot(CfinalNorm,-z_, 'DisplayName',...
%     join(['D = ',num2str(D*1e6),'µm ; rhop = ',num2str(rhoP),'kg.m⁻³']))
% legend('Location','best')
% ylim([-L 0])
% % xlim([0 10])
% grid on
% inv = z(end:-1:1);
% yticks(-inv)
% xlabel("Normed concentration (mps.m⁻³)")
% ylabel("Depth (m)")
% ttl = join(['D = ',num2str(D*1e6),'µm ; rhop = ',num2str(rhoP),'kg.m⁻³']);
% title(ttl)
% hold off

% figure(2)
% plot(dK,-z_)
% 
% figure(3)
% plot(K,-z_)

% % path = '../../Ian/Results/varMP/';
% prof_name = ['profile_D',num2str(D), '_rhop',num2str(rhoP),...
%     '_nPart',num2str(nPart), '_dt', num2str(dt),...
%     '_tf', num2str(tf), '_dtest', num2str(dt_test),'_wind', num2str(WindSpeed),...
%     '_month', num2str(month)];
% exportgraphics(prof,[path,prof_name,'.eps'],'ContentType','vector');
% savefig(prof,[path,'fig/',prof_name,'.fig']);

end

