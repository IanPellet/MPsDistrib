% 
% function [T, MSEh] = MP_transport(dt, nPart, tf)
% Lagrangian transportation model 
% dt : time step
% nPart : total number of particles in the model
% tf : time of simulation
% cKs : coefficient multipling Ks

global L
clear Concentration err

%% Particules initialisation
nPart = 100; % number of particles
zPart = linspace(0, L, nPart); % depth of particles
sizePart = linspace(250e-6, 5e-3, nPart); % size of particles
typePart = fix(1 + (4-1).*rand(nPart)); % type of particles

% allocate memory to store particles
clear part
part(nPart) = MP(typePart(nPart),sizePart(nPart), zPart(nPart), row); % array of MP objects
% Fill the array
for i = 1:(nPart-1)
    part(i) = MP(typePart(i),sizePart(i), zPart(i), row);
end

% Compute initial concentrations
h_init = histogram([part.z_], "BinEdges", z, 'Visible', 'off').Values;
C_init = h_init/dz;


%% Simulation
dt = 1; % Time step
tf = 60*60*12; % Simulation time
nStep = tf/dt+1; % Number of steps
if cast(nStep, 'uint32') ~= nStep
    nStep = cast(nStep+1, 'uint32');
end
T = ones(nStep,1); % allocate memory to store time
C_hist = ones(nStep,N); % allocate memory to store concentration

z_past = [part.z_];
step = 1;
C_hist(step,:) = [C_init];

t=0; OnContinue=true;
T(step,1) = t;
while OnContinue
    step = step+1;
    % Time update
    t=t+dt;
    T(step,1) = t;
    
    %% Particules update
    index = max(1, fix([part.z_]/dz));
    Upart = [part.U_];
    [part.z_] = Step_Lagrangien([part.z_], Upart(index), [part.z_], [part.z_], dt);
    index = max(1, cast(part(1,:)/dz, 'uint32'));
    part(2,:) = u(index);
    part(3,:) = K(index);
    part(4,:) = dK(index);
    
    
    if true

        z_present = part(1,:);
        
        % Computation of the concentration of MPs in each mesh
        h_past = histogram(z_past, "BinEdges", z, 'Visible', 'off').Values;
        h_present = histogram(z_present, "BinEdges", z, 'Visible', 'off').Values;
        CpastNorm = h_past/dz*N/nPart;
        CpresentNorm = h_present/dz*N/nPart;
        C_hist(step,:) = CpresentNorm;
        
        MSEh(step,1) = MSE(CpresentNorm,test);
        
                
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

