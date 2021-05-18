% function [T, MSEh] = MP_transport(dt, nPart, tf)

clear



%% Water column parameters
L = 50;
N = 50;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

%% Diffusivity
DensiteFevrierRhoma
KZ_day = KZ_Fev10;
Row_day = Row_Fev10;
z_day = z_Fev10;
z__day = z__Fev10;
addpath('./fig-Visser/');
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day);

rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 

%% Particules initialisation
nPart = 5e3; % number of particles
zPart = linspace(0, L, nPart); % depth of particles
sizePart = linspace(250e-6, 5e-3, nPart); % size of particle
rhop = 1000;

% allocate memory to store particles
clear part
mp(nPart) = MP; % array of MP objects
% Fill the array
for i = 1:(nPart)
    mp(i) = MP(sizePart(i), rhop, rhow);
end

%% Simulation
dt = 10; % Time step
tf = 60*60; % Simulation time
nStep = tf/dt+1; % Number of steps
if cast(nStep, 'uint32') ~= nStep
    nStep = cast(nStep+1, 'uint32');
end
T = ones(nStep,1); % allocate memory to store time

z_past = zPart;
step = 1;

t=0; OnContinue=true;
T(step,1) = t;
while OnContinue
    step = step+1;
    % Time update
    t=t+dt;
    T(step,1) = t;
    
    %% Particules update
    index = max(1, cast(zPart/dz, 'uint32'));
    U = [mp.U_];
    uz = diag(U(index,:))';
    zPart = Step_Lagrangien(zPart, uz, K(index), dK(index), dt, L);
    
    
    if true

        z_present = zPart;
        
        % Computation of the concentration of MPs in each mesh
        h_past = histogram(z_past, "BinEdges", z, 'Visible', 'off').Values;
        h_present = histogram(z_present, "BinEdges", z, 'Visible', 'off').Values;
        CpastNorm = h_past/dz*N/nPart;
        CpresentNorm = h_present/dz*N/nPart;

                
        dC = max(abs(CpresentNorm - CpastNorm)/dt);
        
        if (t>=tf)
           OnContinue = false;
        end
        
        if mod(t,60)<dt/2
            disp([' Temps : ' num2str(t/tf*100) '% - Ecart : ' num2str(dC)])
        end
        
        z_past = z_present;
    end
end

