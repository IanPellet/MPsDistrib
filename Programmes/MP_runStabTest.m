clear
runID = strrep(join(split(num2str(clock)),'-'),'.','_'); % unique ID to name files
runID = runID{:};

%% Set parameters
dt_test = 60*60*2;
date = datetime(1999,12,19);
tf = 60*60*24*3;

% Multinet sample characteristics 
Station = 'RN2';
% Date = datetime('3/18/2021');

%% Water column parameters
% Find depth at station RN2
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
% Load file with indexes corresponding to each stations
stationFile = '../Data/stationIJ_CEREGE.mat';
load(stationFile,'stationIJ');
% Get the indexes of the right station
I0 = stationIJ{stationIJ{:,'station'} == Station,'I0'};
J0 = stationIJ{stationIJ{:,'station'} == Station,'J0'};
% Get depth at RN2
load(SauvegardeModeleHydro, 'H0')
L = H0(I0,J0);
clear H0,

N = 50;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

nPart = 10e3;
sizeP = 350e-6;
zPart = sparse(linspace(0, L, nPart)); % depth of particles

wind = 10;
frag = 0;

% wind_test = [0 50 100];
% rhop_test = 1000:10:1050;
rhop_test = 1010;

%% Run simulations
nRes = length(rhop_test);
Resultats(nRes) = struct('rhop', NaN, 'zFinal', NaN(1,N), 'dt', NaN);
for iRhop = 1:length(rhop_test)
    rhop = rhop_test(iRhop);
    
    % Diffusivity
    [KZ_day,Row_day,z_day,z__day] = KsSalTemp(wind, date);
    [K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day');
    rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 
    
    % create list of particles
    mp = getMPlist(nPart, ones(nPart,1)*sizeP, rhop, rhow, frag);

    % run simulation
    [saveInt,dt,filePath] = MP_simStabTest(mp, zPart, K, dK, L, dz, tf, dt_test, runID);
end, clear iRhop,

disp('Compute dC')
% compute dC
dtStab = 30*60;
dCinter = 0:dtStab:tf;
dCinter_ = dCinter(1:end-1)+dtStab/2;
timeLine = 0:dt:tf;
dCiBound = NaN(size(dCinter));
for i = 1:length(dCiBound)
    dCiBound(i) = find(timeLine==dCinter(i));
end, clear i,
dCiBound(end) = length(tf);

nCase = dtStab/dt;

dC = NaN(length(dCinter)-2,1);
i = 0;
for iload = 2:length(saveInt)
    file = load([filePath num2str(saveInt(iload)), '.mat'], "zHistory");
    hist = [file.zHistory];
    if isempty(hist{end})
        hist = hist(1:end-1);
    end
    iStart = 0;
    iEnd = 0;
    
    for i30 = 1:length(hist)/nCase
        iStart = iEnd + 1;
        iEnd = iEnd + nCase;
        i = i+1;
        dC(i) = testStability(hist(iStart:iEnd), z, z_, dz, L); 
    end
    clear hist i30,
end, clear i iload,

f1 = figure(1);
plot(dCinter_/60/60,dC*100)
xlabel('simulation time (h)')
ylabel('concentration profile variation (%)')
title('Evolution of the stability of the concentration profile over time',...
    ['Time average on every ' num2str(dtStab/60) ' min of simulation'])


path = '../Results/MP_runStabTest/';

logFile = fopen([path runID '.txt'],'w');
fprintf(logFile, ['Run ID : ' runID '\n'...
    'dt_test = ' num2str(dt_test) ' s \n'...
    'date : ' datestr(date) '\n'...
    'Station : ' Station '\n'...
    'L = ' num2str(L) ' m \n'...
    'N = ' num2str(N) '\n'...
    'nPart = ' num2str(nPart) '\n'...
    'sizeP = ' num2str(sizeP*1e6) ' µm \n'...
    'initial particle repartition : uniform \n'...
    'wind = ' num2str(wind) ' km.h⁻¹ \n'...
    'dtStab = ' num2str(dtStab) ' s']);    
fclose(logFile);
    
figName = [path runID '-dC.fig'];
savefig(f1, figName);



