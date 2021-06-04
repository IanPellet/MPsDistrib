% DensiteFevrierRhoma
% clearvars -except KZ_Fev10 Row_Fev10 z_Fev10 z__Fev10,

dt_test = 60*60*2;
date = datetime(2020,03,18);
% date = "10Fev";
tf = 60*60*24*10;

%% Water column parameters
%% Load hydrodinamic model data
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro, 'H0', 'Lon', 'Lat')

%% Water column parameters
% Find depth of the column
% Lon0 = 5.29; Lat0 = 43.24; %point Somlit
% [I0,J0] = ReperePoint(Lon,Lat,Lon0,Lat0); % data indices corresponding to the location
Station = 'RN2';
% Load file with indexes corresponding to each stations
stationFile = '../Data/stationIJ_CEREGE.mat';
load(stationFile,'stationIJ');
% Get the indexes of the right station
I0 = stationIJ{stationIJ{:,'station'} == Station,'I0'};
J0 = stationIJ{stationIJ{:,'station'} == Station,'J0'};
L = H0(I0,J0); % depth
clear H0 Lon Lat stationIJ,

N = 50;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

nPart = 10e3;
% pd = makedist('Normal', 'mu', 350e-6, 'sigma', 50e-6);
% sizeP = random(pd,nPart,1);
% clear pd,
sizeP = ones(nPart,1)*350e-6;
zPart = linspace(0,L,nPart);
frag = 0;

wind_test = [100];
rhop_test = 1025;

testStab = 60*60*5; % stability test interval
dtAvgC = 30*60;
tinter = 0:dtAvgC:tf;
tC = tinter(1:end-1)+dtAvgC/2;

path = '../Results/MP_runStabTest/';



%% Run simulations
simIDs = cell(length(rhop_test)*length(wind_test),1);
iID = 0;
for iRhop = 1:length(rhop_test)
    rhop = rhop_test(iRhop);
    for iWind = 1:length(wind_test)
        wind = wind_test(iWind);
        
        runID = strrep(join(split(num2str(clock)),'-'),'.','_'); % unique ID to name files
        runID = runID{:};
        iID = iID+1;
        simIDs{iID} = runID;

        % Diffusivity
        [KZ_day,Row_day,z_day,z__day] = KsSalTemp(wind, date);
        
%         KZ_day = KZ_Fev10;
%         Row_day = Row_Fev10;
%         z_day = z_Fev10;
%         z__day = z__Fev10;
      
        [K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day');

        rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 
        clear KZ_day Row_day z_day z__day,
    
        % create list of particles
        mp = getMPlist(nPart, sizeP, rhop, rhow, frag);

%% Run simulation
        [dt, historyFiles] = MP_simStabTest(mp, zPart, K, dK, L, dz, tf, dt_test, runID, 60*60*24);
        
%% Compute concentration profiles time average
        disp('Compute C')
        nCase = dtAvgC/dt;
        meanConc = cell(ceil(tf/dtAvgC),1);
        stdConc = cell(ceil(tf/dtAvgC),1);
        iConc = 0;
        for iFile = 1:length(historyFiles)
            load(historyFiles{iFile}, "zHistory");
            for iHist = 1:length(zHistory)
                if mod(iHist,nCase) == 0
                    iConc = iConc +1;
                    [meanConc{iConc},stdConc{iConc}] = getMeanConc(zHistory(iHist-nCase+1:iHist), N, dz);
                end
            end, clear iHist,
        end, clear iFile iConc zHistory,
        
%% Evaluate stability of the concentration profile
        disp('Compute StabC')
        [StabC,tStabC] = getStability(meanConc, testStab, dtAvgC, tC);
        
%% Save parameters and results to file
        disp('Save parameters and results')
        save([path runID '.mat'],...
            'runID', 'dt_test', 'date', 'tf', 'L', 'N', 'dz', 'nPart',...
            'sizeP', 'zPart', 'frag', 'wind', 'rhop', 'testStab', 'path',...
            'historyFiles', 'meanConc', 'stdConc', 'StabC', 'K', 'dK', 'dt',...
            'testStab', 'tC', 'tStabC', 'dtAvgC');
        
%% Plot figures
        f1 = figure(1);
        plot(tStabC/60/60,StabC*100)
        xlabel('simulation time (h)')
        ylabel('Root Mean Square Error (mps.m⁻³)')
        title('Evolution of the stability of the concentration profile over time',...
            ['Tested between profiles separated by ' num2str(testStab/60/60) ' h of simulation'])
        
        meanC = cell2mat(meanConc);
        f2 = figure(2); clf,
        h = pcolor(tC/60/60',-z_,meanC');
        a = colorbar;
        a.Label.String = 'Concentration (mps.m⁻³)';
        xlabel('Simulation time (h)')
        ylabel('Depth (m)')
        title('Evolution of the concentration profile over time',...
            ['Time average on every ' num2str(dtAvgC/60) ' min of simulation'])
        set(h, 'EdgeColor', 'none');
        
        f3 = figure(3); clf,
        plot(K,-z_)
        xlabel('Diffusivity (m².s⁻¹)')
        ylabel('Depth (m)')
        title('Diffusivity profile');

%% Save figures
        figName = [path runID '-dC.fig'];
        savefig(f1, figName);
        figName = [path runID '-C.fig'];
        savefig(f2, figName);
        figName = [path runID '-K.fig'];
        savefig(f3, figName);
        
        clear dC figName historyFiles meanC meanConc
        
    end, clear iWind,
end, clear iRhop iID,

% hold on
% for iID = 1:length(simIDs)
%     dC = load([path 'dC/' simIDs{iID} '-dC.mat'], 'dC');
%     plot(dtC/60/60,dC.dC*100, 'DisplayName', simIDs{iID})
% end, clear iID,
% legend('Location', 'best')
% xlabel('simulation time (h)')
% ylabel('concentration profile variation (%)')
% title('Evolution of the stability of the concentration profile over time',...
%     ['Time average on every ' num2str(dtStab/60) ' min of simulation'])
% hold off
% 
% runID = strrep(join(split(num2str(clock)),'-'),'.','_'); % unique ID to name files
% runID = runID{:};
% fig2Name = [path runID '-dCcomparison.fig'];
% savefig(f2, fig2Name);