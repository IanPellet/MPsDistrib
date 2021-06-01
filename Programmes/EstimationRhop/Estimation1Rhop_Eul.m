function [Resultats,ConcentrationSample,Ztest_,z_] = Estimation1Rhop_Eul(SizePart, RhoP_test, saveFig, N, windSpeed)
%%ESTIMATION1RHOP_EUL Finds the modeled particle density closest
%%to data
% SizePart : (double) size of the modeled particles (m)
% RhoP_tes : (double array) modeled particle densities tested (kg.m⁻³)
% saveFig : (boolean) if true figures will be saved to '../Results/EstimRho/'
% N : (int) default = 500, number of meshes. High : +precise -fast / Low : -precise +fast
% windSpeed : (double) wind speed (km.h⁻¹)
% 
% example : Estimation1Rhop_Eul(400e-6, 900:1100, true, 600)

%% Set Parameters
path = '../Results/EstimRho/'; % saved figures directory

% Multinet sample characteristics 
Station = 'RN2';
Date = datetime('3/18/2021');

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
% L = 50;

if nargin < 4
    N = 500; % number of meshes. High : +precise -fast / Low : -precise +fast
end

dz= L/N; % size of meshes
z = (0:dz:L)'; % meshes boundaries
z_ = z(1:end-1)+dz/2; % center of the meshes

% day = '10fev'; % day corresponding to diffusive turbulence data
day = false;

%% Data

dh = 0.71; % Net oppening (m)

% Concentration at each depth and filtered volume
ConcVolFile = '../Data/ConcVol_MP.txt';
ConcVolTable = load_ConcVol_data(ConcVolFile);

% Multinet condition
concCond = strcmp(ConcVolTable{:,'station'},Station) & ConcVolTable{:,'date'} == Date;

% Concentration data
CMes = ConcVolTable{concCond,'C'};
ZMes = ConcVolTable{concCond,'depth'}+dh/2;

% CMes=[0.62 0.34 0.06 0.02 0]'; % Concentration
% ZMes=[1 10 15 40 50]'+dh/2; % Depth of the samples
ConcentrationSample = CMes;
DepthSample = ZMes;
% DataInterp = interp1(ZMes,CMes,z_,'pchip')'; % Interpolated data


%% Find corresponding depths to get from the model
boundTest = zeros(length(DepthSample)*2,1); % Boundaries in between which modeled part number have to be tested
ibound = zeros(size(boundTest)); % index boundaries of meshes to avg

for i = 1:length(boundTest)
    if mod(i,2)
        boundTest(i) = DepthSample(fix((i+1)/2))-dh/2;
    else
        boundTest(i) = DepthSample(fix((i+1)/2))+dh/2;
    end
    ibound(i) = min(length(z_),max(1,fix(boundTest(i)/dz)+1));
end, clear i,

Ztest_ = ZMes; 

%% Simulations
clear Resultats
% Initalisation of results data structure fields
Resultats(1:length(RhoP_test)) = struct('RhoP', 0, 'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'rmseErreur', 0);

if SizePart
    modSize = SizePart;
else
    modSize = 350e-6;
end

% launch simulation for each rhop
for iRes = 1:length(RhoP_test)
    RhoP = RhoP_test(iRes);    
    
    [conc, z_] = Transport_Eulerian(modSize, RhoP, N, L, day, windSpeed, Date);
%     [conc, z_] = Transport_Eulerian(modSize, RhoP, N, L, day);
    
    % find the modeled concentration at depth corresponding with sample
    ConcentrationModel = zeros(size(ConcentrationSample));
    j = 0;
    for ic = 1:2:length(ibound)
        j = j+1;
        ConcentrationModel(j) = mean(conc(ibound(ic):ibound(ic+1)));
    end, clear ic j,
    
    % Compute the proportionnality coefficient minimizing the square error
    % between model and data
    alpha = ConcentrationModel\ConcentrationSample;
    
    % Compute the root mean square error between model and data
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
    
    % Fill the result structure 
    Resultats(iRes).RhoP = RhoP;
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).Erreur = Erreur;
    Resultats(iRes).rmseErreur = rmseErreur;
    Resultats(iRes).conc = conc';
end, clear iRes,

% Find rhop correponding to minimal error
minI = [Resultats.rmseErreur] == min([Resultats.rmseErreur]); 
minRho = Resultats(minI).RhoP;

%% Display results
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm']; % Figures titles

% figure 1 : plot all simulation results
f1 = figure(1); clf,
% plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
xlim([0 3])
ylim([-L 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
for res = Resultats
    plot(res.conc*res.Alpha,-z_,'DisplayName', ['RhoP = ' num2str(res.RhoP) 'kg.m⁻³'])
end
legend('Location', 'best')
title(ttl)
hold off

% figure 2 : plot error with respect to rhoP
f2 = figure(2); clf,
hold on
plotErr = zeros(size([Resultats]));
for i = 1:length(Resultats)
    plotErr(i) = Resultats(i).rmseErreur;
end
plot([Resultats.RhoP], plotErr)
xlabel('Rho_p (kg.m⁻³)');
ylabel('RMSE (mps.m⁻³)');
title(ttl)
hold off

% figure 3 : plot best rhoP profile
f3 = figure(3); clf,
% plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
xlim([0 3])
ylim([-L 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on 
plot(CMes,-Ztest_,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
plot(Resultats(minI).conc*Resultats(minI).Alpha,-z_,'DisplayName', ['RhoP = ' num2str(Resultats(minI).RhoP) 'kg.m⁻³'])
legend('Location', 'best')
title([ttl ' -- Rho_p = ' num2str(minRho) 'kg.m⁻³'])
hold off

%% Save figures
if saveFig
    if length(RhoP_test) == 1
        rhoInter = num2str(RhoP_test);
    else
        rhoInter = [num2str(min(RhoP_test)) '-' num2str(RhoP_test(2)-RhoP_test(1)) '-' num2str(max(RhoP_test))];
    end
    
    fileName = ['size' num2str(modSize*1e6) '_rho' rhoInter '_dataMarine'];
    
    F = [f1, f2, f3];
    xPart = '1part_';
    N = {'profils_', 'error_', 'min_'};
    for i=1:length(F)
        f = F(i);
        n = N{i};
        NamePNG = [path 'png/' xPart n fileName '.png'];
        NameFIG = [path 'fig/' xPart n fileName '.fig'];
        exportgraphics(f,NamePNG);
        savefig(f,NameFIG);
    end, clear i,
end

end
