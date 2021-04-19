function [T, errorh] = Ajust_model(dt, nPart, tf, Ks, v, dt_test)
% Lagrangian transportation model 
% dt : time step
% nPart : total number of particles in the model
% tf : time of simulation
% Ks : coefficient multipling Ks

global L
clear Concentration err

%% Speed computation formulas
Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];
indNom=3;

%% Load hydrodinamic model
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)

%% Water column parameters
L = 50;
N=50;  dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh
z0=L*Sigma;   

%% Particules initialisation
D=350e-6; % Particle diametre
rop=1010.5; % Particle density

z_part = linspace(0, L, nPart);

%% Speed initialisation
DensiteFevrierRhoma 
% [K,dK] = wcp_interpolation(z0,KZ_Fev10*cKs,-z_); % Diffusivity
K_val = Ks;
K = ones(size(z))*K_val;
dK = zeros(size(z));


InitialisationVitesseTransport
S=rop./row;     
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
%u=Ws; u(rop<row)=-Ws(rop<row);
% u(1) = 0; u(end) = 0;
u = -ones(size(Ws))*v;

%% Particules matrix definition
% part = [z1 z2 ... zn ; 
%         u1 u2 ... un ;
%         K1 K2 ... Kn ;
%        dK1 dK2...dK3]
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
errorh = ones(nStep,1);
test = C_analytical(v, Ks, z_, nPart, L)*N/nPart;

z_past = part(1,:);
step = 1;
ChNorm(step,:) = [CiNorm];
errorh(step,1) = totalError(CiNorm,test);

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
        
        errorh(step,1) = totalError(CpresentNorm,test);
        
                
        dC = max(abs(CpresentNorm - CpastNorm)/dt);
        
        if (t>=tf)
           OnContinue = false;
        end
        
        if mod(t,60)<dt/2
            disp([' Temps : ' num2str(t/tf*100) '% - Ecart : ' num2str(dC/sum(CpresentNorm)*100), '%'])
        end
        
        z_past = z_present;
    end
end

errorh = errorh/N;

%% Plot end profile
prof = figure(1);
clf
% hold on

% ax1 = subplot(1,2,1);
% hold on
% plot(ax1, CiNorm,-z_, 'DisplayName','t = 0')
% plot(ax1, CpresentNorm,-z_, 'DisplayName',['t = ', num2str(tf)])
hold on
plot(test,-z_, 'DisplayName','Analytical Solution')
plot(CpresentNorm,-z_, 'DisplayName','Lagrangian Model')
legend('Location','best')
ylim([-L 0])
xlim([0 5])
grid on
inv = z(end:-1:1);
yticks(-inv)
xlabel("Normed concentration (mps.m⁻³)")
ylabel("Depth (m)")
ttl = join(["Ks = ", num2str(Ks), "m².s⁻¹"]);
title(ttl)

% ax2 = subplot(1,2,2);
% plot(ax2, K(1:end-1), -z, 'DisplayName',['Ks = ', num2str(Ks)])
% ylim([-L 0])
% grid on
% inv = z(end:-1:1);
% yticks(-inv)
% xlabel("Turbulent diffusion (m².s⁻¹)")
% legend('Location','best')


hold off
prof_name = ['../../Ian/Results/dt', num2str(dt), '_nPart', num2str(nPart),...
    '_tf', num2str(tf), '_Ks', num2str(Ks), '_v', num2str(v),'_dtest', num2str(dt_test), '_profile.eps'];
exportgraphics(prof,prof_name,'ContentType','vector');
