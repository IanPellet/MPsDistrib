function [zFinal, dt, partList, tauxAgg, runID, AggContent] = Part_Simulator2(partList, partZinit, K, dK, H, dz, rhof, tf, dt_test, dtTheo, S, Buoy, Agg, PlotPart, saveHist,saveLastSec,addAggT,runID,nmp)
%PART_SIMULATOR2 Run simulation with particles
%
% PARAMETERS
% partList : MP or aggregate objects array, modeled particles
% partZinit : double array, initial particle's position (m)
% K : double array, diffusivity profile (m².s⁻¹)
% dK : double array, diffusivity gradient (m.s⁻¹)
% H : double, water column depth (m)
% dz : double, water column discretisation step (m)
% rhof : fluid (sea water) density 
% tf : double, simulation time (s)
% dt_test : double, time interval at witch to display simulation advencement (s)
% dtTheo : time step initialisation (might be changed if too big) (seconds)
% S :
% Buoy : boolean, simulate particle buoyancy (true) or not (false) ?
% Agg : boolean, simulate particle aggregation (true) or not (false) ?
% PlotPart : boolean, plot particle position during simulation (true) or 
%            not (false) ?
% saveHist : boolean, save particle position/composition during simulation 
%            (true) or not (false) ?
% saveLastSec : double, period over which to save data (s), if you want to 
%               save the data over the whole simulation time set it to tf, 
%               if you want to only save the last minute of simulation (for 
%               example) set it to 60             
% addAggT : double, time at which to launch aggregation (seconds) you can 
%           start a simulation without particle aggregation, then start 
%           aggregation at time addAggT
% runID : str, ID given to this simulation 
% nmp : boolean, save the number of MP contained in each aggregate at each
%       time step
%
% OUTPUT
% zFinal : double array, final particle's position (m)
% dt : double, time step used for simulation (s)
% partList : MP or aggregate objects array, modeled particles
% tauxAgg :
% runID : str, ID given to this simulation 
% AggContent :


saveTo = '../Results/AggPlot/'; % where to save data

% saveLastSec = tf; 

fprintf(['\n\n--------------------- Simulation running ---------------------\n'])
saveTime = 60*60*24;
savei = 1;

% default arguments
if nargin < 16
    nmp = false;
    saveLastSec = tf;
    addAggT = 0;
end
if nargin < 15
    saveHist = true;
end
if nargin < 12
    Buoy = true;
    Agg = true;
    PlotPart = true;
end

K = reshape(K,1,[]);
dK = reshape(dK,1,[]);
    
    partZ = reshape(partZinit, 1, numel(partZinit)); % Particles positions (m)
    clear partZinit,
    FreePart = [partList.Locked]==0;
    
%     d = 1;
    StepPD = makedist('Normal', 'mu', 0, 'sigma', 1);
% StepPD = makedist('Uniform', 'lower', -1, 'upper', 1);
    d = std(StepPD)^2;
    
    if Buoy
        partI= cast(max(1, fix(partZ/dz)), 'uint32'); % int array, index of each particle's current mesh
        partWs = VitesseNguyen(rhof, partList, partI); % MP fall velocities on the column (m.s⁻¹)
        
%         if ~Agg & 
            Aire = [partList.A];
            Peri = [partList.P];
            Vol = [partList.V];
            rho = [partList.Rho];
%             L = [partList.L];
            partWs = VitesseNguyen(rhof, partList, partI, Aire, Peri, Vol, rho);
%             partWs = wsAhrens(rhof, partList, partI, L, rho);
%         end
    else
        partWs = zeros(size(partList));
    end
    
    nMP = sum([partList.Fb]==0);
    P = length(partList);
    
    %% Time step initialisation
    if nargin >= 11
        dt = dtTheo;
    else
        dt = 60; % double, time step (s)
    end
    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
   
    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs(partWs)))); % check condition dt < dz/max|u|
    dt = min(min(dt));
    
    clear ddK,

    %% Init history 
%     saveHist = nargin > 8 && saveLastSec ~= 0;
    if saveHist
        timeToSave = min(saveLastSec,saveTime);
        saveNstep = ceil(timeToSave/dt)+1;
        zHistory = cell(saveNstep,1);
        tauxAgg = nan(saveNstep,1);
        
        saveStep = 1;
        zHistory{saveStep} = partZ;
        
        
        nMPAgg = P-sum(FreePart);
        tauxAgg(saveStep) = nMPAgg/nMP;
    end
    if nmp
        AggContent = cell(saveNstep,1);
        AggContent{saveStep} = [partList.nContent];
    end
    
    %% Simulation
    t=0; 
    OnContinue=true;
        
    while OnContinue
        % Plot
        if PlotPart
            plotPartState(partList, partZ, H, t)
        end
        
        % Aggregation
        if Agg && t>=addAggT
            [partList, FreePart] = Part_Aggregation(partList, FreePart, partZ, dK, K, dt, dz, S, Buoy, rhof);
        end
        
        %% Time update
        t=t+dt;
        
        % Particules update
        partI= cast(max(1, fix(partZ/dz)), 'uint8'); % int array, index of each particle's current mesh
        
        if Buoy
            % Find current fall velocity of each particle
            if Agg && t>=addAggT
                partWs = VitesseNguyen(rhof, partList, partI); % MP fall velocities on the column (m.s⁻¹)
            else
                partWs = VitesseNguyen(rhof, partList, partI, Aire, Peri, Vol, rho);
%                 partWs = wsAhrens(rhof, partList, partI, L, rho);
            end
        end
            
        % Update position of aggregates and free MPs
        partZ(FreePart) = Step_Lagrangien_noDist(partZ(FreePart), partWs(FreePart), K(partI(FreePart)), dK(partI(FreePart)), dt, H, StepPD);
        
        % Update position of locked MPs (same position as the aggregate
        % it's locked on)
        partZ(~FreePart) = partZ([partList(~FreePart).Locked]);
        

      
        %%
        % Save to history
        if saveHist && t >= tf-saveLastSec
            saveStep = saveStep+1;
            zHistory{saveStep} = partZ;
            nMPAgg = P-sum(FreePart);
            tauxAgg(saveStep) = nMPAgg/nMP;
        end
        if nmp
            AggContent{saveStep} = [partList.nContent];
        end
        if saveTime*savei<=t && t<tf
            if nmp
                save([saveTo runID{:} '_' num2str(savei) '.mat'], 'zHistory', 'tauxAgg', 'AggContent')
            else
                save([saveTo runID{:} '_' num2str(savei) '.mat'], 'zHistory', 'tauxAgg')
            end
            savei = savei+1;
            
            timeToSave = min(saveLastSec-t,saveTime);
            saveNstep = ceil(timeToSave/dt)+1;
            zHistory = cell(saveNstep,1);
            tauxAgg = nan(saveNstep,1);
        
            saveStep = 0;
        end
    
        % Test
%         if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 )
        if rem(t,dt_test)+dt > dt_test     
%             % Display progress
%             if mod(t,60)<dt/2
                disp([' Temps : ' num2str(t/tf*100) '%'])
%             end
        end
        % Test if final time is reached
        if (t>=tf)
           OnContinue = false;
        end

    end
    if saveHist
        zFinal = zHistory; % 2D double array, final particle's position (m)
    else
        zFinal = {partZ}; % 1D double array, final particle's position (m)
    end
    if nmp
        save([saveTo runID{:} '_' num2str(savei) '.mat'], 'zHistory', 'tauxAgg', 'AggContent')
    else
        save([saveTo runID{:} '_' num2str(savei) '.mat'], 'zHistory', 'tauxAgg')
    end
    
    
    % Plot
    if PlotPart
        plotPartState(partList, partZ, H, t)
    end

end

function plotPartState(partList, partZ, H, t)
    nPart = length(partZ);
    scatter(1:nPart, -partZ, [partList.L]*1e5/2, [partList.Fb],'filled')
    xlim([1 nPart])
    ylim([-H 0])
    xlabel("Particles")
    ylabel("Depth (m)")
    hh = fix(t/60/60);
    mm = fix(t/60-hh*60);
    title(['t = ', num2str(hh), ':', num2str(mm)])
    pause(0)
end



