function [dt, historyFiles] = MP_simStabTest(mp, zInit, K, dK, L, dz, tf, dt_test, ID, dtSave)
%MP_SIMULATOR Run simulation with particles in mp
%
% PARAMETERS
% mp : MP objects array, modeled particles
% zInit : double array, initial particle's position (m)
% K : double array, diffusivity profile (m².s⁻¹)
% dK : double array, diffusivity gradient (m.s⁻¹)
% L : double, water column depth (m)
% dz : double, water column discretisation step (m)
% tf : double, simulation time (s)
% dt_test : double, test interval (s)
% saveLastSec : double, time over which to save particle's positions
% history (s) //Optional

% dtSave : int, saving time intevals (s)
%
% OUTPUT
% dt : int, simulation's time step (s)
% historyFiles : string cell array, files containing concentration history
%

path = ['../Results/MP_runStabTest/Data/' ID '-zHist'];
clear ID,

fprintf(['\n\n--------------------- nPart = ' num2str(length(mp)) ' ---------------------\n'])

%% Init particle speed and position
U = [mp.U_]; % 2D double array, mp fall velocities on the column (m.s⁻¹)   
uz = NaN(size(mp)); % double array, mp fall velocities (m.s⁻¹)
zPart = reshape(zInit,1,numel(zInit)); % double array, particle's position (m)
clear zInit mp,


%% Time step initialisation
dt = 10; % double, time step (s)
ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
dt = min(dt, dz/max(max(abs(U)))); % check condition dt < dz/max|u|
clear ddK,


%% Init history
zHistory = cell(ceil(dtSave/dt),1);
historyFiles = cell(ceil(tf/dtSave),1);

saveTime = dtSave;
saveStep = 0;
iFile = 0;

%% Simulation
t=0; 
OnContinue=true;
while OnContinue
    % Time update
    t=t+dt;

    %% Particules update
    index = cast(max(1, fix(zPart/dz)), 'uint8'); % int array, index of each particle's current mesh
    % Find current fall velocity of each particle
    for i=1:length(zPart)
        uz(i) = U(index(i),i);
    end
    % Update particle's postion
    zPart = Step_Lagrangien(zPart, uz, K(index), dK(index), dt, L);

    % Save to history
    saveStep = saveStep+1;
    zHistory{saveStep} = zPart;
    
    if t >= tf || abs(t-saveTime) <= dt/2
        disp("Saving history...");

        saveTime = saveTime + dtSave; 
        iFile = iFile +1;
        saveStep = 0;
        
        historyFiles{iFile} = [path num2str(t), '.mat'];
        save(historyFiles{iFile}, "zHistory");

    end

    % Test
    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 )
        % Test if final time is reached
        if (t>=tf)
           OnContinue = false;
        end
        % Display progress
        if mod(t,60)<dt/2
            disp([' Temps : ' num2str(t/tf*100) '%'])
        end
    end

end
    
end