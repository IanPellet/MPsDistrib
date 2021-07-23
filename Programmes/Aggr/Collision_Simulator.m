function [collision, dt] = Collision_Simulator(mpList, aggList, mpZinit, aggZinit, K, dK, L, dz, tf, dt_test, dtTheo)


fprintf(['\n\n--------------------- Simulation running ---------------------\n'])

    figure(5); clf,

    mpU = [mpList.Ws]; % MP fall velocities on the column (m.s⁻¹) 
    aggU = [aggList.Ws]; % OrgaAggr fall velocities on the column (m.s⁻¹) 
    
    mpZ = reshape(mpZinit, 1, numel(mpZinit)); % MP positions (m)
    aggZ = reshape(aggZinit, 1, numel(aggZinit)); % OrgaAggr positions (m)
    clear mpZinit aggZinit,
    
    mpUz = NaN(size(mpList)); % mp fall velocities at depth z (m.s⁻¹)
    aggUz = NaN(size(aggList)); % OrgaAggr fall velocities at depth z (m.s⁻¹)

    sizePart = [mpList.Size aggList.Size];
    nPart = length(sizePart);

    SumD = (reshape([mpList.Size], 1, length(mpList)) + reshape([aggList.Size], length(aggList), 1));
    
    col = [repmat([0 0.4470 0.7410],length(mpZ),1) ; repmat([0.4660 0.6740 0.1880],length(aggZ),1)];
    
    collision = cell(0);
    
    %% Time step initialisation

    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
    
    ampu = mat2cell(mpU, [length(mpList(1).Ws)], ones(1,length(mpList)));
    baggu = mat2cell(aggU, [length(aggList(1).Ws)], ones(1,length(aggList)));
    B = 6.*sqrt(2.*K');
    D = mpList(1).Size + aggList(1).Size;

    dtmax = nan(length(ampu), length(baggu));
    for ia = 1:length(ampu)
        for ib = 1:length(baggu)
            cdelu = abs(ampu{ia}-baggu{ib});
            cdelu_ = (cdelu(1:end-1)+cdelu(2:end))./2;
            det = 72.*K' + 4.*cdelu_.*D;
            X = [(-B + sqrt(det))./(2.*cdelu_) (-B - sqrt(det))./(2.*cdelu_)];
            X = X(X>=0);
            dtmax(ia,ib) = min(X.^2);
            

        end, clear ib,
    end, clear ia,
    
    dt = min(min(dtmax));
    clear ampu baggu cdelu dtmax X det,

    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs([mpU aggU])))); % check condition dt < dz/max|u|
    
    clear ddK meanDU,

    %% Simulation
    t=0; 
    boom = 0; % number of collisions 
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
        
        mpZn_1 = mpZ;
        aggZn_1 = aggZ;
        % Update position 
        mpZ = Step_Lagrangien(mpZ, mpUz, K(mpI), dK(mpI), dt, L);
        aggZ = Step_Lagrangien(aggZ, aggUz, K(aggI), dK(aggI), dt, L);
        
        %% Collision
        mpPos = reshape(mpZ, 1, length(mpZ));
        aggPos = reshape(aggZ, length(aggZ), 1);
        DeltaPos = abs(mpPos - aggPos);
      
        
        [aggClose, mpClose] = find(DeltaPos <= SumD./2);
         for iClose = 1:length(aggClose)
%              disp('POUET')
            iagg = aggClose(iClose);
            imp = mpClose(iClose);
            
            boom = boom +1;
            
            mpI= cast(max(1, fix(mpZn_1(imp)/dz)), 'uint8'); % int array, index of MP's mesh at t-dt
            aggI= cast(max(1, fix(aggZn_1(iagg)/dz)), 'uint8'); % int array, index of Agg's mesh at t-dt
            % Find fall velocity of each particle
            mpUz2 = mpU(mpI, imp);
            aggUz2 = aggU(aggI, iagg);
            
            delU = mpUz2 - aggUz2;
            deldK = dK(mpI)-dK(aggI)';
            sumSqrtK = sqrt(K(mpI)) + sqrt(K(aggI)');
            delMax = abs(delU + deldK).*dt + sumSqrtK.*sqrt(6*dt);
            
            collision{boom} = [mpZn_1(imp)-aggZn_1(iagg) delMax (mpList(imp).Size+aggList(iagg).Size)/2];
            
        end


        

%         %% Plot 
%         zPart = [mpZ aggZ];
%         scatter(1:length(zPart), -zPart, sizePart*1e5/2, col,'filled')
%         xlim([1 nPart])
%         ylim([-31 -30])
%         xlabel("Particles")
%         ylabel("Depth (m)")
% %         hh = fix(t/60/60);
% %         mm = fix(t/60-hh*60);
% %         title(['t = ', num2str(hh), ':', num2str(mm)])
%         title(['t = ', num2str(t)])
%         pause(0)

    
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