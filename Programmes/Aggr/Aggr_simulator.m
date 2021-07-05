function [zFinal, dt] = Aggr_simulator(mp, xInit, yInit, zInit, K, dK, L, dz, tf, dt_test, saveLastSec)
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

fprintf(['\n\n--------------------- nPart = ' num2str(length(mp)) ' -- save last ' num2str(saveLastSec) 's ---------------------\n'])

    U = [mp.U_]; % 2D double array, mp fall velocities on the column (m.s⁻¹)   
    uz = NaN(size(mp)); % double array, mp fall velocities (m.s⁻¹)
    
    xPart = reshape(xInit,1,numel(xInit)); % double array, particle's position (m)
    yPart = reshape(yInit,1,numel(yInit)); % double array, particle's position (m)
    zPart = reshape(zInit,1,numel(zInit)); % double array, particle's position (m)
    clear xInit yInit zInit,
    saveHist = nargin > 8 && saveLastSec ~= 0;
    
    %% Time step initialisation
    dt = 10; % double, time step (s)
    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs(U)))); % check condition dt < dz/max|u|
    clear ddK,
    
    %% Init fragmentation rate list
    rFragMP = [mp.fragRate_]; % Tester si on a eu fragmentation et récup que si c'est le cas
    if unique(rFragMP) == 0
        frag = false;
        clear rFragMP,
    else
        frag = true;
    end
    if frag
        rFragDT = rFragMP*dt; % Get frangmentation rate for the time dt
    end
    %% Init history 
    if saveHist
        saveNstep = ceil(saveLastSec/dt);
        History = cell(saveNstep,3);
        saveStep = 0;
    end
    
    %% Simulation
    t=0; 
    OnContinue=true;
    while OnContinue
        % Time update
        t=t+dt;
        
        %% Fragment particles
        if frag
            draw = rand(size(rFragMP));
            mpFrag = rFragDT >= draw;
            if sum(mpFrag) > 0
                [mp, zPart] = MP_fragmentation(mp, zPart, mpFrag);
                U = [mp.U_]; % 2D double array, mp fall velocities on the column (m.s⁻¹)
                rFragMP = [mp.fragRate_]; % Tester si on a eu fragmentation et récup que si c'est le cas
                rFragDT = rFragMP*dt; % Get frangmentation rate for the time dt
            end
        end
        %% Particules update
        index = cast(max(1, fix(zPart/dz)), 'uint8'); % int array, index of each particle's current mesh
        % Find current fall velocity of each particle
        for i=1:length(zPart)
            uz(i) = U(index(i),i);
        end
        % Update particle's postion
        [xPart, yPart, zPart] = Step_Aggr(xPart, yPart, zPart, uz, K(index), dK(index), dt, L);
        
        dX = diff(sort(xPart));
        dY = diff(sort(yPart));
        dZ = diff(sort(zPart));
        
        if sum(dX<500e-6 & dY<500e-6 & dZ<500e-6) ~= 0
%             disp('POUET')
        end
        
        % Save to history
        if saveHist && t >= tf-saveLastSec
            saveStep = saveStep+1;
            History{saveStep,1} = xPart;
            History{saveStep,2} = yPart;
            History{saveStep,3} = zPart;
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
    if saveHist
        zFinal = History; % 2D double array, final particle's position (m)
    else
        zFinal = [xPart;yPart;zPart]; % 1D double array, final particle's position (m)
    end
    
end