function [zFinal, dt] = Part_Simulator(mpList, aggList, mpZinit, aggZinit, K, dK, L, dz, tf, dt_test)


fprintf(['\n\n--------------------- Simulation running ---------------------\n'])

    mpU = [mpList.Ws]; % MP fall velocities on the column (m.s⁻¹) 
    aggU = [aggList.Ws]; % OrgaAggr fall velocities on the column (m.s⁻¹) 
    
    mpZ = reshape(mpZinit, 1, numel(mpZinit)); % MP positions (m)
    aggZ = reshape(aggZinit, 1, numel(aggZinit)); % OrgaAggr positions (m)
    clear mpZinit aggZinit,
    
    mpUz = NaN(size(mpList)); % mp fall velocities at depth z (m.s⁻¹)
    aggUz = NaN(size(aggList)); % OrgaAggr fall velocities at depth z (m.s⁻¹)
    
%     part = [mpList aggList];
    sizePart = [mpList.Size aggList.Size];
    nPart = length(sizePart);
    
    %% Time step initialisation
    dt = 60; % double, time step (s)
    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs([mpU aggU])))); % check condition dt < dz/max|u|
    clear ddK,

%     %% Init history 
%     saveHist = nargin > 8 && saveLastSec ~= 0;
%     if saveHist
%         saveNstep = ceil(saveLastSec/dt);
%         zHistory = cell(saveNstep,1);
%         saveStep = 0;
%     end
    
    %% Simulation
    t=0; 
    OnContinue=true;
    while OnContinue
        % Time update
        t=t+dt;
        
        %% Particules update
        mpI= cast(max(1, fix(mpZ/dz)), 'uint8'); % int array, index of each particle's current mesh
        aggI= cast(max(1, fix(aggZ/dz)), 'uint8'); % int array, index of each particle's current mesh
        
        % Find current fall velocity of each particle
        for i=1:length(mpZ)
            mpUz(i) = mpU(mpI(i),i);
        end
        for i=1:length(aggZ)
            aggUz(i) = aggU(aggI(i),i);
        end
        
        % Update particle's postion
        mpUnLock = [mpList.Locked] == 0;
        
        mpZ(mpUnLock) = Step_Lagrangien(mpZ(mpUnLock), mpUz(mpUnLock), K(mpI(mpUnLock)), dK(mpI(mpUnLock)), dt, L);
        aggZ = Step_Lagrangien(aggZ, aggUz, K(aggI), dK(aggI), dt, L);
        
        %% Aggregation
        zPart = [mpZ aggZ];
        
        % Sort particles by position
        [zPartSort, iSort] = sort(zPart);
        rayonSort = sizePart(iSort)/2;
%         partSort = part(iSort);
        
        % Compute particles distances
        movSumRayon2 = movsum(rayonSort,2);
        minEquart = movSumRayon2(2:end);
        diffZpart = diff(zPartSort);
        
        iClose = find(diffZpart<=minEquart);
        for i=iClose
            % 1 AGG PEUT INTÉGRER PLUSIEURS MP À CHAQUE STEP PAS ENCORE IMPLEMENTÉ 
            % ADHERENCEN VARIABLE DES AGGRÉGATS NON IMPLÉMENTÉE
            if (iSort(i)<=length(mpList) && iSort(i+1)>length(mpList))...
                    || (iSort(i+1)<=length(mpList) && iSort(i)>length(mpList))
%                 partSort(i).Type ~= partSort(i+1).Type 
                iagg = max(iSort(i),iSort(i+1))-length(mpList); % find the index of the aggregate
                imp = min(iSort(i),iSort(i+1)); % find the index of the MP
                
                [aggList(iagg), mpList(imp)] = aggList(iagg).aggrMP(mpList(imp));
            end
        end
        
        %% Plot
        scatter(1:length(zPart), -zPart, sizePart*1e4)
        xlim([1 nPart])
        ylim([-L 0])
        xlabel("Particles")
        ylabel("Depth (m)")
        pause(0)
        
%         % Save to history
%         if saveHist && t >= tf-saveLastSec
%             saveStep = saveStep+1;
%             zHistory{saveStep} = zPart;
%         end
    
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
%     if saveHist
%         zFinal = zHistory; % 2D double array, final particle's position (m)
%     else
%         zFinal = zPart; % 1D double array, final particle's position (m)
%     end
    zFinal = {mpZ ; aggZ};
end