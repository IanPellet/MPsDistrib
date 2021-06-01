function [saveInter, dt, path] = MP_simStabTest(mp, zInit, K, dK, L, dz, tf, dt_test, ID)
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
%
% OUTPUT
% zFinal : double array, final particle's position (m)
%

path = ['../Results/MP_runStabTest/Data/' ID '-zHist'];

fprintf(['\n\n--------------------- nPart = ' num2str(length(mp)) ' ---------------------\n'])

    U = [mp.U_]; % 2D double array, mp fall velocities on the column (m.s⁻¹)   
    uz = NaN(size(mp)); % double array, mp fall velocities (m.s⁻¹)
    zPart = zInit; % double array, particle's position (m)
    clear zInit,

    
    %% Time step initialisation
    dt = 10; % double, time step (s)
    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs(U)))); % check condition dt < dz/max|u|
    clear ddK,
    

    %% dC parm
    dtSave = 60*60*24;
    saveInter = 0:dtSave:tf;
%     dCinter_ = dCinter(1:end-1)+dtStab/2;
    timeLine = 0:dt:tf;
    iBoundSave = NaN(size(saveInter));
    for i = 1:length(iBoundSave)
        iBoundSave(i) = find(timeLine==saveInter(i));
    end, clear i,
    iBoundSave(end) = tf;

    
    %% Init history 
    zHistory = cell(iBoundSave(2),1);
    saveStep = 0;
    iTest = 2;

    %% Simulation
    t=0; 
    OnContinue=true;
    while OnContinue
        % Time update
        t=t+dt;
        
    
        %% Particules update
        index = max(1, fix(zPart/dz)); % int array, index of each particle's current mesh
        % Find current fall velocity of each particle
        for i=1:length(zPart)
            uz = U(index(i),i);
        end
        % Update particle's postion
        zPart = Step_Lagrangien(zPart, uz, K(index), dK(index), dt, L);
        
        % Save to history
        saveStep = saveStep+1;
        zHistory{saveStep} = zPart;
        
        if t >= tf || abs(t-saveInter(iTest)) <= dt/2
            disp("Saving history...");
            iTest = iTest+1;
            saveStep = 0;
            save([path num2str(t), '.mat'], "zHistory");
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

    clearvars -except saveInter dt path
    
end