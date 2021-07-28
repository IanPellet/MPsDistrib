function [collision, dt,mpZ,aggZ,mpList,aggList] = Collision_Simulator(mpList, aggList, mpZinit, aggZinit, K, dK, L, dz, tf, dt_test, dt_coll)


fprintf(['\n\n--------------------- Simulation running ---------------------\n'])

    figure(1); clf,

    mpU = [mpList.Ws]; % MP fall velocities on the column (m.s⁻¹) 
    aggU = [aggList.Ws]; % OrgaAggr fall velocities on the column (m.s⁻¹) 
%     mpU = mpU.*0;
    aggU = aggU.*0;
    
    mpZ = reshape(mpZinit, 1, numel(mpZinit)); % MP positions (m)
    aggZ = reshape(aggZinit, 1, numel(aggZinit)); % OrgaAggr positions (m)
    clear mpZinit aggZinit,
    
    mpUz = NaN(size(mpList)); % mp fall velocities at depth z (m.s⁻¹)
    aggUz = NaN(size(aggList)); % OrgaAggr fall velocities at depth z (m.s⁻¹)
%     mpUz = zeros(size(mpUz));
    aggUz = zeros(size(aggUz));

    sizePart = [mpList.Size aggList.Size];
    nPart = length(sizePart);
    
    mpSize = reshape([mpList.Size], 1, length(mpList));
    aggSize = reshape([aggList.Size], length(aggList), 1);
    
    mpFree = [mpList.Locked]==0;
    mpIndex = 1:length(mpList);
    
    col = [repmat([0 0.4470 0.7410],length(mpZ),1) ; repmat([0.4660 0.6740 0.1880],length(aggZ),1)];
    
    collision = cell(0);
    
    
    d = 1;
    StepPD = makedist('Normal', 'mu', 0, 'sigma', d);
    
    %% Time step initialisation

    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
    Rmax = d*4;
    D = mpList(1).Size + aggList(1).Size;
    
%     ampu = mat2cell(mpU, [length(mpList(1).Ws)], ones(1,length(mpList)));
%     baggu = mat2cell(aggU, [length(aggList(1).Ws)], ones(1,length(aggList)));
%     
%     B = 6.*Rmax.*sqrt(2/d.*K');
%     
% 
%     dtmax = nan(length(ampu), length(baggu));
%     for ia = 1:length(ampu)
%         for ib = 1:length(baggu)
%             cdelu = abs(ampu{ia}-baggu{ib});
%             cdelu_ = (cdelu(1:end-1)+cdelu(2:end))./2;
%             
%             det = 4.*Rmax^2.*2/d.*K' + 4.*cdelu_.*D;
%             X = [(-B + sqrt(det))./(2.*cdelu_) (-B - sqrt(det))./(2.*cdelu_)];
%             X = X(X>=0);
%             dtmax(ia,ib) = min(X.^2);
%             
% 
%         end, clear ib,
%     end, clear ia,
%     
%     dt = min(min(dtmax));
%     clear ampu baggu cdelu dtmax X det,

    dt = min((D./(2.*Rmax.*sqrt(2/d.*K'))).^2);

    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs([mpU aggU])))); % check condition dt < dz/max|u|
    
    clear ddK meanDU B D,
    
    %% INIT collisions
    boom = 0; % number of collisions 
%     X = round(dt_coll/dt);
%     sauveXdt = cell(X,2);
%     t_Xdti = -X; % Index at which we can find the positions at t-Xdt
%     iSauve = 0; % Index to at which to save curent values

    %% Simulation
    t=0;
    X = 0;
    if t+dt/2 >= X*dt_coll
%         disp('POUET')
        X = X+1;
        mpt_Xdt = mpZ;
        aggt_Xdt = aggZ;
    end
    
    
    
    
    %% Collision
        mpPos = reshape(mpZ(mpFree), 1, length(mpZ(mpFree)));
        aggPos = reshape(aggZ, length(aggZ), 1);
        DeltaPos = abs(mpPos - aggPos);
      
        SumD = (mpSize(mpFree) + aggSize);
        
        [aggClose, mpClose] = find(DeltaPos <= SumD./2);
        
        mpIndexFree = mpIndex(mpFree);
        
%         if ~ isempty(aggClose) && t_Xdti>0
        if ~ isempty(aggClose)
            % get mp and agg positions at t-Xdt
%             mpt_Xdt = sauveXdt{t_Xdti,1};
%             aggt_Xdt = sauveXdt{t_Xdti,2};
            
            for iClose = 1:length(aggClose)
                iagg = aggClose(iClose);
                imp = mpIndexFree(mpClose(iClose));
                
                boom = boom +1;

                collision{boom} = [mpt_Xdt(imp)-aggt_Xdt(iagg) (mpList(imp).Size+aggList(iagg).Size)/2];
                
                [aggList(iagg), mpList(imp)] = aggList(iagg).aggrMP(mpList(imp));
                mpFree(imp) = false; 
            end
        end
    
%     t_Xdti = (t_Xdti>0)*(mod(t_Xdti, X)+1) + (t_Xdti<=0)*(t_Xdti+1);
%     iSauve = mod(iSauve, X) +1;
%     sauveXdt{iSauve,1} = mpZ;
%     sauveXdt{iSauve,2} = aggZ;
           
    
    OnContinue=true;
%       %% Plot
% 
%         zPart = [mpZ aggZ];
%         scatter(1:length(zPart), -zPart, sizePart*1e5/2, col,'filled')
%         xlim([1 nPart])
%         ylim([-31 -30])
%         xlabel("Particles")
%         ylabel("Depth (m)")
%         hh = fix(t/60/60);
%         mm = fix(t/60-hh*60);
%         title(['t = ', num2str(hh), ':', num2str(mm)])
%         pause(0)

        
    while OnContinue
        % Time update
        t=t+dt;
        
%         t_Xdti = (t_Xdti>0)*(mod(t_Xdti, X)+1) + (t_Xdti<=0)*(t_Xdti+1);
%         iSauve = mod(iSauve, X) +1;
    
        
        %% Particules update
        mpI= cast(max(1, fix(mpZ/dz)), 'uint8'); % int array, index of each particle's current mesh
        aggI= cast(max(1, fix(aggZ/dz)), 'uint8'); % int array, index of each particle's current mesh
        
        % Find current fall velocity of each particle
        for i=1:length(mpZ)
            mpUz(i) = mpU(mpI(i),i);
        end
%         for i=1:length(aggZ)
%             aggUz(i) = aggU(aggI(i),i);
%         end
        
        % Update position of aggregates and free MPs 
        mpZ(mpFree) = Step_Lagrangien_noDist(mpZ(mpFree), mpUz(mpFree), K(mpI(mpFree)), dK(mpI(mpFree)), dt, L, StepPD);
%         aggZ = Step_Lagrangien_noDist(aggZ, aggUz, K(aggI), dK(aggI), dt, L, StepPD);
        
        % Update position of locked MPs (same position as the aggregate
        % it's locked on)
        mpZ(~mpFree) = aggZ([mpList(~mpFree).Locked]);
        
        % Save position
%         sauveXdt{iSauve,1} = mpZ;
%         sauveXdt{iSauve,2} = aggZ;
        if t+dt/2 >= X*dt_coll
%             disp('POUET')
            X = X+1;
            mpt_Xdt = mpZ;
            aggt_Xdt = aggZ;
        end

        %% Collision
        mpPos = reshape(mpZ(mpFree), 1, length(mpZ(mpFree)));
        aggPos = reshape(aggZ, length(aggZ), 1);
        DeltaPos = abs(mpPos - aggPos);
      
        SumD = (mpSize(mpFree) + aggSize);
        
        [aggClose, mpClose] = find(DeltaPos <= SumD./2);
        
        mpIndexFree = mpIndex(mpFree);
        
%         if ~ isempty(aggClose) && t_Xdti>0
        if ~ isempty(aggClose)
            % get mp and agg positions at t-Xdt
%             mpt_Xdt = sauveXdt{t_Xdti,1};
%             aggt_Xdt = sauveXdt{t_Xdti,2};
            
            for iClose = 1:length(aggClose)
                iagg = aggClose(iClose);
                imp = mpIndexFree(mpClose(iClose));
                
                boom = boom +1;

                collision{boom} = [mpt_Xdt(imp)-aggt_Xdt(iagg) (mpList(imp).Size+aggList(iagg).Size)/2];
                
                [aggList(iagg), mpList(imp)] = aggList(iagg).aggrMP(mpList(imp));
                mpFree(imp) = false; 
            end
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
    
     %% Plot
        zPartPlot = [mpZ aggZ];
        
        if sum(~mpFree)~=0
            col(~mpFree,:) = repmat([0.8500 0.3250 0.0980], sum(~mpFree),1);
            col(length(mpZ) + [mpList(~mpFree).Locked], :) = repmat([0.4940 0.1840 0.5560], sum(~mpFree),1);
        end
        
        
        scatter(1:length(zPartPlot), -zPartPlot, sizePart*1e5/2, col,'filled')
        xlim([1 nPart])
%         ylim([-max(zPartPlot) -min(zPartPlot)])
        xlabel("Particles")
        ylabel("Depth (m)")
%         hh = fix(t/60/60);
%         mm = fix(t/60-hh*60);
%         title(['t = ', num2str(hh), ':', num2str(mm)])
        title(['Real collisions -- t = ', num2str(t)])
        pause(0)
    

end