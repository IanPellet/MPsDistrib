function [zFinal, dt, mpList, aggList] = Part_Simulator(mpList, aggList, mpZinit, aggZinit, K, dK, L, dz, tf, dt_test, dtTheo)


fprintf(['\n\n--------------------- Simulation running ---------------------\n'])
    
    saveHist = true;
    saveLastSec = tf;

    path = '../Results/Aggr/Video/';
    fileName = 'aggrGauss';
    f1 = figure(5); clf,

    mpU = [mpList.Ws]; % MP fall velocities on the column (m.s⁻¹) 
    aggU = [aggList.Ws]; % OrgaAggr fall velocities on the column (m.s⁻¹) 
    
    mpZ = reshape(mpZinit, 1, numel(mpZinit)); % MP positions (m)
    aggZ = reshape(aggZinit, 1, numel(aggZinit)); % OrgaAggr positions (m)
    clear mpZinit aggZinit,
    
    mpUz = NaN(size(mpList)); % mp fall velocities at depth z (m.s⁻¹)
    aggUz = NaN(size(aggList)); % OrgaAggr fall velocities at depth z (m.s⁻¹)

    sizePart = [mpList.Size aggList.Size];
    nPart = length(sizePart);
    
    mpSize = reshape([mpList.Size], 1, length(mpList));
    aggSize = reshape([aggList.Size], length(aggList), 1);
    
    mpFree = [mpList.Locked]==0;
    mpIndex = 1:length(mpList);
    
    col = [repmat([0 0.4470 0.7410],length(mpZ),1) ; repmat([0.4660 0.6740 0.1880],length(aggZ),1)];
    
    %% Time step initialisation
    if nargin == 11
        dt = dtTheo;
    else
        dt = 60; % double, time step (s)
    end
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

    %% Init history 
%     saveHist = nargin > 8 && saveLastSec ~= 0;
    if saveHist
        saveNstep = ceil(saveLastSec/dt);
        zHistory = cell(saveNstep,1);
        saveStep = 0;
    end
    
    %% Simulation
    t=0; 
    OnContinue=true;
      %% Plot
        col = [repmat([0 0.4470 0.7410],length(mpZ),1) ; repmat([0.4660 0.6740 0.1880],length(aggZ),1)];
        if sum(~mpFree)~=0
            col(~mpFree,:) = repmat([0.8500 0.3250 0.0980], sum(~mpFree),1);
            col(length(mpZ) + [mpList(~mpFree).Locked], :) = repmat([0.4940 0.1840 0.5560], sum(~mpFree),1);
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
%         exportgraphics(f1, [path fileName num2str(t) '.png']);
        
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
        
        
        % Update position of aggregates and free MPs
        mpZ(mpFree) = Step_Lagrangien(mpZ(mpFree), mpUz(mpFree), K(mpI(mpFree)), dK(mpI(mpFree)), dt, L);
        aggZ = Step_Lagrangien(aggZ, aggUz, K(aggI), dK(aggI), dt, L);
        
        % Update position of locked MPs (same position as the aggregate
        % it's locked on)
        mpZ(~mpFree) = aggZ([mpList(~mpFree).Locked]);
        
%         %% Aggregation
%         closeParticles = true;
%         
%         while closeParticles && t>=60+dt
%             mpFree = [mpList.Locked] == 0; % isfree MPs
%             
%             zPart = [mpZ(mpFree) aggZ]; % get position of aggregates and free MPs
% 
%             % Sort particles by position
%             [zPartSort, iSort] = sort(zPart);
%             rayonSort = sizePart(iSort)/2;
% 
%             % Compute particles distances
%             movSumRayon2 = movsum(rayonSort,2);
%             minEquart = movSumRayon2(2:end);
%             diffZpart = diff(zPartSort);
% 
%             iClose = find(diffZpart<=minEquart);
%             for i=iClose
%                 % 1 AGG PEUT INTÉGRER PLUSIEURS MP À CHAQUE STEP PAS ENCORE IMPLEMENTÉ 
%                 % ADHERENCEN VARIABLE DES AGGRÉGATS NON IMPLÉMENTÉE
%                 iagg = max(iSort(i),iSort(i+1))-length(mpList); % find the index of the aggregate
%                 imp = min(iSort(i),iSort(i+1)); % find the index of the MP
%                 if iagg>0 && imp<=length(mpList) && mpList(imp).Locked ==0
%                     test = abs(aggZ(iagg)-mpZ(imp)) <= (aggList(iagg).Size+mpList(imp).Size)/2;
%                     disp(test)
%                     if ~test
%                         error('Pouet')
%                     end
%                     [aggList(iagg), mpList(imp)] = aggList(iagg).aggrMP(mpList(imp));
% %                     disp(['AGG ' num2str(iagg) ' ' num2str(imp)])
%                 else
%                     closeParticles = false;
%                 end
%             end
%             if isempty(iClose)
%                 closeParticles = false;
%             end
%         end
        %% Aggregation
        mpPos = reshape(mpZ(mpFree), 1, length(mpZ(mpFree)));
        aggPos = reshape(aggZ, length(aggZ), 1);
        DeltaPos = abs(mpPos - aggPos);

        mpSizeFree = mpSize(mpFree);
%         SumR = (mpSizeFree + aggSize)./2;
        SumD = (mpSizeFree + aggSize);
        
        mpI= cast(max(1, fix(mpPos/dz)), 'uint8'); % int array, index of each MP's current mesh
        aggI= cast(max(1, fix(aggPos/dz)), 'uint8'); % int array, index of each Agg's current mesh
        % Find current fall velocity of each particle
        mpUz2 = nan(size(mpPos));
        for i=1:length(mpPos)
            mpUz2(i) = mpU(mpI(i),i);
        end
        
        aggUz2 = nan(size(aggPos));
        for i=1:length(aggPos)
            aggUz2(i) = aggU(aggI(i),i);
        end
        delU = mpUz2 - aggUz2;
        deldK = dK(mpI)-dK(aggI)';
        sumSqrtK = sqrt(K(mpI)) + sqrt(K(aggI)');
        delMax = abs(delU + deldK).*dt + sumSqrtK.*sqrt(6*dt);
%         DelMean = abs(delU + deldK).*dt;
        
%         a = -1./DelMax;
%         b = 1 + SumR./DelMax;
%         
%         p = a.*DeltaPos+b;

        c2 = 0.5;
        d = 1./(1-exp(delMax.*(delMax+SumD)./(2*c2)));
        a = (1-d).*exp(SumD.^2/(8*c2));
        p = a.*exp(-DeltaPos.^2./(2*c2))+d;
        

%         [aggClose, mpClose] = find(DeltaPos <= SumR);
        [aggClose, mpClose] = find(p > 0);
        
        mpIndexFree = mpIndex(mpFree);

        for iClose = 1:length(aggClose)
            iagg = aggClose(iClose);
            imp = mpIndexFree(mpClose(iClose));
            
%             if mpList(imp).Locked ~= 0
%                 error('Locked')
%             end
%             if abs(mpZ(imp)-aggZ(iagg)) > (mpList(imp).Size+aggList(iagg).Size)/2
%                 error('Far')
%             end
            
            if p(aggClose(iClose),mpClose(iClose))>=rand
                [aggList(iagg), mpList(imp)] = aggList(iagg).aggrMP(mpList(imp));
                mpFree(imp) = false;
            end
        end

        %% Plot
        zPartPlot = [mpZ aggZ];
        
        if sum(~mpFree)~=0
            col(~mpFree,:) = repmat([0.8500 0.3250 0.0980], sum(~mpFree),1);
            col(length(mpZ) + [mpList(~mpFree).Locked], :) = repmat([0.4940 0.1840 0.5560], sum(~mpFree),1);
        end
        
        
        scatter(1:length(zPartPlot), -zPartPlot, sizePart*1e5/2, col,'filled')
        xlim([1 nPart])
        ylim([-L 0])
        xlabel("Particles")
        ylabel("Depth (m)")
        hh = fix(t/60/60);
        mm = fix(t/60-hh*60);
        title(['t = ', num2str(hh), ':', num2str(mm)])
        pause(0)
%         exportgraphics(f1, [path fileName num2str(t) '.png']);
        %%
        % Save to history
        if saveHist && t >= tf-saveLastSec
            saveStep = saveStep+1;
            zHistory{saveStep} = zPartPlot;
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
        zFinal = zHistory; % 2D double array, final particle's position (m)
    else
        zFinal = {mpZ ; aggZ}; % 1D double array, final particle's position (m)
    end
    

end