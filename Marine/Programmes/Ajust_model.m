function [T, MSEh] = Ajust_model(dt, nPart, tf)
% Lagrangian transportation model 
% 
% 

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

z_part = z_' * ones(1,nPart);
z_part = reshape(z_part, size(z_part,2)*size(z_part,1),1);
z_part = sort(z_part)';

%% Speed initialisation
DensiteFevrierRhoma 
[K,dK] = wcp_interpolation(z0,KZ_Fev10,-z_); % Diffusivity
%K_val = 0.01;
%K = ones(size(z))*K_val;
%dK = zeros(size(z));
InitialisationVitesseTransport
S=rop./row;     
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
%u=Ws; u(rop<row)=-Ws(rop<row);
u = Ws*0;

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
disp(class(h_init))
CiNorm = h_init/dz/nPart;

if false
    figure(1), clf
    ax1 = subplot(1,3,1);
    plot(ax1,-part(1,:),'x')
    ylim([-L 0])
    grid on
    inv = z(end:-1:1);
    yticks(-inv)

    lim = [0.5 1.5];
    ax2 = subplot(1,3,2);
    plot(ax2,CiNorm,-z_)
    xlim(lim)
    ylim([-L 0])
    grid on
    yticks(-inv)

    ax3 = subplot(1,3,3);
    plot(ax3,(2/d*K(1:end-1)*dt).^(1/2), -z)
    ylim([-L 0])
    grid on
    yticks(-inv)
    pause(1)
end

%% Simulation
nStep = tf/dt+1;
if cast(nStep, 'uint32') ~= nStep
    nStep = cast(nStep+1, 'uint32');
end
T = ones(nStep,1);
ChNorm = ones(nStep,N);
MSEh = ones(nStep,1);
test = ones(size(CiNorm));

z_past = part(1,:);
step = 1;
ChNorm(step,:) = [CiNorm];
MSEh(step,1) = MSE(CiNorm,test);

t=0; OnContinue=true;
T(step,1) = t;
while OnContinue
    step = step+1;
    % Time update
    t=t+dt;
    T(step,1) = t;
    % Particules update
    part(1,:) = Step_Lagrangien(part(1,:), part(2,:), part(3,:), part(4,:), dt);
    index = max(1, cast(part(1,:)/dz, 'uint32'));
    part(2,:) = u(index);
    part(3,:) = K(index);
    part(4,:) = dK(index);
    
    if true

        z_present = part(1,:);
        
        % Computation of the concentration of MPs in each mesh
        h_past = histogram(z_past, "BinEdges", z, 'Visible', 'off').Values;
        h_present = histogram(z_present, "BinEdges", z, 'Visible', 'off').Values;
        CpastNorm = h_past/dz/nPart;
        CpresentNorm = h_present/dz/nPart;
        ChNorm(step,:) = CpresentNorm;
        
        MSEh(step,1) = MSE(CpresentNorm,test);
        
                
        dC = max(abs(CpresentNorm - CpastNorm)/dt);
        
        if (t>=tf)
           OnContinue = false;
        end
        
        if mod(t,30)<dt/2
            disp([' Temps : ' num2str(t/tf*100) '% - Ecart : ' num2str(dC)])
        end
        
        if false
            ax1 = subplot(1,3,1);
            plot(ax1,-part(1,:),'x')
            ylim([-L 0])
            grid on
            inv = z(end:-1:1);
            yticks(-inv)

            ax2 = subplot(1,3,2);
            plot(ax2,CpresentNorm,-z_)
            ylim([-L 0])
            xlim(lim)
            grid on
            inv = z(end:-1:1);
            yticks(-inv)

            ax3 = subplot(1,3,3);
            plot(ax3,(2/d*K(1:end-1)*dt).^(1/2), -z)
            ylim([-L 0])
            grid on
            inv = z(end:-1:1);
            yticks(-inv)

            pause(0.1)
        end
        
        z_past = z_present;
    end
end
%error = 
%Delta = [T mean(ChNorm,2) min(ChNorm,[],2) max(ChNorm,[],2)];
