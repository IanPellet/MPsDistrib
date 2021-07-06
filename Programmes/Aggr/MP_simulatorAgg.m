function [zFinal, dt, mp] = MP_simulatorAgg(mp, zInit, K, dK, L, dz, tf, dt_test, saveLastSec)
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

    zPart = reshape(zInit,1,numel(zInit)); % double array, particle's position (m)
    clear zInit,
    saveHist = nargin > 8 && saveLastSec ~= 0;
    
    %% Aggr
    sizePart = [mp.size_];
    nPart = length(mp);
    
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
        zHistory = cell(saveNstep,1);
        saveStep = 0;
    end
    
    %% Simulation
    t=0; 
    OnContinue=true;
    figure(1), clf,
    
    
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
        uz = NaN(size(mp)); % double array, mp fall velocities (m.s⁻¹)
        for i=1:length(zPart)
            uz(i) = U(index(i),i);
        end
        % Update particle's postion
        zPart = Step_Lagrangien(zPart, uz, K(index), dK(index), dt, L);
        
        %% Aggregation
        [zPartSort, iSort] = sort(zPart);
        rayonSort = sizePart(iSort)/2;
        movSumRayon2 = movsum(rayonSort,2);
        minEquart = movSumRayon2(2:end);
        diffZpart = diff(zPartSort);
        
        iClose = find(diffZpart<=minEquart);
        iAgg = zeros(size(mp));
        for i=iClose
            if rand < mp(i+1).sticky_
                mp(i+1).size_ = mp(i).size_ + mp(i+1).size_;
                iAgg(i) = 1;
            end
        end
        mp = mp(~iAgg);
        zPart = zPart(~iAgg);
        
        % Save to history
        if saveHist && t >= tf-saveLastSec
            saveStep = saveStep+1;
            zHistory{saveStep} = zPart;
        end
        
        scatter(1:length(zPart), -zPart, [mp.size_]*1e4)
        xlim([1 nPart])
        ylim([-L 0])
        xlabel("Particles")
        ylabel("Depth (m)")
        pause(0)

    
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
        zFinal = zHistory; % 2D double array, final particle's position (m)
    else
        zFinal = zPart; % 1D double array, final particle's position (m)
    end
    
end