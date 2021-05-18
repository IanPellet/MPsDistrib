% function [T, MSEh] = MP_transport(dt, nPart, tf)

clear

tf = 60*60*24*2; % Simulation time
dt_test = 60*60;

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
sizePart = linspace(250e-6, 500e-6, nPart); % size of particle
% sizePart = ones(size(zPart))*400e-6;
rhop = 1020;

% allocate memory to store particles
clear part
mp(nPart) = MP; % array of MP objects
% Fill the array
for i = 1:(nPart)
    mp(i) = MP(sizePart(i), rhop, rhow);
end

U = [mp.U_];

%% Time step init
dt = 10;
ddK = diff(dK)/dz;
dt = min(dt, abs(min(1./ddK)/10)); % condition dt<<min(1/ddK) 
dt = min(dt, dz/max(max(abs(U)))); % condition dt < dz/max|u|

%% Initial conditions
h_init = histogram(zPart, "BinEdges", z, 'Visible', 'off').Values;
CiNorm = h_init/dz*N/nPart;

%% Simulation
nStep = tf/dt_test+1;
if cast(nStep, 'uint32') ~= nStep
    nStep = cast(nStep+1, 'uint32');
end
T = ones(nStep,1);
ChNorm = ones(nStep,N);

z_past = zPart;
step = 1;
ChNorm(step,:) = [CiNorm];

t=0; OnContinue=true;
T(step,1) = t;
% uz = zeros(size(zPart));
while OnContinue
    % Time update
    t=t+dt;
    
    %% Particules update
    index = max(1, cast(zPart/dz, 'uint32'));

    for i=1:length(zPart)
        uz = U(index(i),i);
    end
    
    
    zPart = Step_Lagrangien(zPart, uz, K(index), dK(index), dt, L);
    

    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 ) 
    %if true
        step = step+1;
        T(step,1) = t;
        z_present = zPart;
        
        % Computation of the concentration of MPs in each mesh
        h_past = histogram(z_past, "BinEdges", z, 'Visible', 'off').Values;
        h_present = histogram(z_present, "BinEdges", z, 'Visible', 'off').Values;
        CpastNorm = h_past/dz*N/nPart;
        CpresentNorm = h_present/dz*N/nPart;
        ChNorm(step,:) = CpresentNorm;     
                
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
CfinalNorm = CpresentNorm;
plot(CfinalNorm, -z_)

