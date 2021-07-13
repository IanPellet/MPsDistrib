function [zFinal, dt, mpList, aggList] = Part_Simulator(mpList, aggList, mpZinit, aggZinit, K, dK, L, dz, tf, dt_test, dtTheo)


fprintf(['\n\n--------------------- Simulation running ---------------------\n'])

    path = '../Results/Aggr/Video/';
    f1 = figure(1); clf,

    mpU = [mpList.Ws]; % MP fall velocities on the column (m.s⁻¹) 
    aggU = [aggList.Ws]; % OrgaAggr fall velocities on the column (m.s⁻¹) 
    
    mpZ = reshape(mpZinit, 1, numel(mpZinit)); % MP positions (m)
    aggZ = reshape(aggZinit, 1, numel(aggZinit)); % OrgaAggr positions (m)
    clear mpZinit aggZinit,
    
    mpUz = NaN(size(mpList)); % mp fall velocities at depth z (m.s⁻¹)
    aggUz = NaN(size(aggList)); % OrgaAggr fall velocities at depth z (m.s⁻¹)

    sizePart = [mpList.Size aggList.Size];
    nPart = length(sizePart);
    
    %% Time step initialisation
    dt = 10; % double, time step (s)
    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
    
%     cumulU = max(abs(mpU-aggU),[], 2);
%     cumulU_ = (cumulU(1:end-1) + cumulU(2:end))/2;
%     cumulK_ = (2.*sqrt(6.*K))';
%     cumulUK_ = cumulU_ + cumulK_;
%     dt = max((min([mpList.Size]) + min([aggList.Size]))./cumulUK_);
    
    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs([mpU aggU])))); % check condition dt < dz/max|u|
    
    clear ddK,
% dt = dtTheo;

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
      %% Plot
        mpUnLock = [mpList.Locked] == 0;
        mpLock = ~mpUnLock;
        col = [repmat([0 0.4470 0.7410],length(mpZ),1) ; repmat([0.4660 0.6740 0.1880],length(aggZ),1)];
        if sum(mpLock)~=0
            col(mpLock,:) = repmat([0.8500 0.3250 0.0980], sum(mpLock),1);
            col(length(mpZ) + [mpList(mpLock).Locked], :) = repmat([0.4940 0.1840 0.5560], sum(mpLock),1);
        end
        zPart = [mpZ aggZ];
        scatter(1:length(zPart), -zPart, sizePart*1e5/2, col,'filled')
        xlim([1 nPart])
        ylim([-L 0])
        xlabel("Particles")
        ylabel("Depth (m)")
        hh = fix(t/60/60);
        mm = fix(t/60-hh*60);
        title(['t = ', num2str(hh), ':', num2str(mm)])
        pause(0)
%         exportgraphics(f1, [path 'aggrdt' num2str(t) '.png']);
        
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
        mpLock = ~mpUnLock;
        
        mpZ(mpUnLock) = Step_Lagrangien(mpZ(mpUnLock), mpUz(mpUnLock), K(mpI(mpUnLock)), dK(mpI(mpUnLock)), dt, L);
        aggZ = Step_Lagrangien(aggZ, aggUz, K(aggI), dK(aggI), dt, L);
        
        mpZ(mpLock) = aggZ([mpList(mpLock).Locked]);
        
        %% Aggregation
        closeParticles = true;
        
        while closeParticles && t>=60+dt
            zPart = [mpZ([mpList.Locked]==0) aggZ];

            % Sort particles by position
            [zPartSort, iSort] = sort(zPart);
            rayonSort = sizePart(iSort)/2;

            % Compute particles distances
            movSumRayon2 = movsum(rayonSort,2);
            minEquart = movSumRayon2(2:end);
            diffZpart = diff(zPartSort);

            iClose = find(diffZpart<=minEquart);
            for i=iClose
                % 1 AGG PEUT INTÉGRER PLUSIEURS MP À CHAQUE STEP PAS ENCORE IMPLEMENTÉ 
                % ADHERENCEN VARIABLE DES AGGRÉGATS NON IMPLÉMENTÉE
                iagg = max(iSort(i),iSort(i+1))-length(mpList); % find the index of the aggregate
                imp = min(iSort(i),iSort(i+1)); % find the index of the MP
                if ((iSort(i)<=length(mpList) && iSort(i+1)>length(mpList))...
                        || (iSort(i+1)<=length(mpList) && iSort(i)>length(mpList)))...
                        && mpList(imp).Locked ==0
                    
                    

                    [aggList(iagg), mpList(imp)] = aggList(iagg).aggrMP(mpList(imp));
%                     disp(['AGG ' num2str(iagg) ' ' num2str(imp)])
                else
                    closeParticles = false;
                end
            end
            if isempty(iClose)
                closeParticles = false;
            end
        end
        
        %% Plot
        mpUnLock = [mpList.Locked] == 0;
        mpLock = ~mpUnLock;
        col = [repmat([0 0.4470 0.7410],length(mpZ),1) ; repmat([0.4660 0.6740 0.1880],length(aggZ),1)];
        if sum(mpLock)~=0
            col(mpLock,:) = repmat([0.8500 0.3250 0.0980], sum(mpLock),1);
            col(length(mpZ) + [mpList(mpLock).Locked], :) = repmat([0.4940 0.1840 0.5560], sum(mpLock),1);
        end
        
        zPartPlot = [mpZ aggZ];
        scatter(1:length(zPartPlot), -zPartPlot, sizePart*1e5/2, col,'filled')
        xlim([1 nPart])
        ylim([-L 0])
        xlabel("Particles")
        ylabel("Depth (m)")
        hh = fix(t/60/60);
        mm = fix(t/60-hh*60);
        title(['t = ', num2str(hh), ':', num2str(mm)])
        pause(0)
%         exportgraphics(f1, [path 'aggrdt' num2str(t) '.png']);
        
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