function [zFinal, dt, mpList, aggList] = Part_Simulator(mpList, aggList, mpZinit, aggZinit, K, dK, L, dz, tf, dt_test, dtTheo)
%PART_SIMULATOR Simulator with aggregation
%
% PARAMETERS
% mpList : MP objects array, modeled MP particles
% aggList : Aggr object array, modeled aggregates
% mpZinit : double array, initial MP positions (m)
% aggZinit : double array, initial Aggr positions (m)
% K : double array, diffusivity profile (m².s⁻¹)
% dK : double array, diffusivity gradient (m.s⁻¹)
% L : double, water column depth (m)
% dz : double, water column discretisation step (m)
% tf : double, simulation time (s)
% dt_test : double, test interval (s)
% dtTheo : double, time step (s) can be changed if needed
%
% OUTPUTS
% zFinal : double array, finial particle position (m)
% dt : double, time step used for simulation (s)
% mpList : MP objects array, modeled MP particles
% aggList : Aggr object array, modeled aggregates

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
        

        %% Aggregation
        
        
        mpPos = reshape(mpZ(mpFree), 1, length(mpZ(mpFree)));
        aggPos = reshape(aggZ, length(aggZ), 1);
        DeltaPos = (mpPos - aggPos);

        SumD = (mpSize(mpFree) + aggSize);
        
        mpI= cast(max(1, fix(mpPos/dz)), 'uint8'); % int array, index of each MP's current mesh
        aggI= cast(max(1, fix(aggPos/dz)), 'uint8'); % int array, index of each Agg's current mesh
        
        % Find current fall velocity of each particle
        mpUz2 = nan(1,sum(mpFree));
        mpi=0;
        for i=1:length(mpFree)
            if mpFree(i)
                mpi = mpi+1;
                mpUz2(mpi) = mpU(mpI(mpi),i);
            end
        end
        
        aggUz2 = nan(size(aggPos));
        for i=1:length(aggPos)
            aggUz2(i) = aggU(aggI(i),i);
        end
        delU = mpUz2 - aggUz2;
        deldK = dK(mpI)-dK(aggI)';

        eqtyp = 1;
        
        delMax = abs(delU + deldK).*dt + 10.*eqtyp.*sqrt(K(mpI)) + sqrt(K(aggI)').*sqrt(6*dt);
        
        
        mu = (delU + deldK).*dt;
        sigma = sqrt(2.*eqtyp.*(K(mpI)+K(aggI)').*dt);
        
        p = nan(size(DeltaPos));
        for i1 = 1:size(DeltaPos,1)
            for i2 = 1:size(DeltaPos,2)
                
                if abs(DeltaPos(i1,i2)) <= 10*abs(delMax(i1,i2)) % only compute proba if the particles are close enought

                    dist_move = makedist('Normal','mu',mu(i1,i2),'sigma',sigma(i1,i2));


                    if DeltaPos(i1,i2) < -SumD(i1,i2)/2

                        p(i1,i2) = 1-cdf(dist_move, -SumD(i1,i2)/2-DeltaPos(i1,i2));
                
                    elseif DeltaPos(i1,i2) > SumD(i1,i2)/2
               
                        p(i1,i2) = cdf(dist_move, SumD(i1,i2)/2-DeltaPos(i1,i2));
                   
                    else
                 
                        p(i1,i2) = 1;
                     
                    end
                end
    
            end
        end
        

        [aggClose, mpClose] = find(p > 0);
        
        mpIndexFree = mpIndex(mpFree);

        for iClose = 1:length(aggClose)
            iagg = aggClose(iClose);
            imp = mpIndexFree(mpClose(iClose));
            
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