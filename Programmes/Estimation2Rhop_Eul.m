function [Resultats,minI,ConcentrationSample,Ztest_,z_] = Estimation2Rhop_Eul(SizePart, RhoP_test, saveFig, N, windSpeed)
%%ESTIMATION2RHOP_EUL Finds the modeled particle density closest
%%to data with two different densities
% SizePart : (double) size of the modeled particles (m)
% RhoP_tes : (double array) modeled particle densities tested (kg.m⁻³)
% saveFig : (boolean) if true figures will be saved to '../Results/EstimRho/'
% N : (int) default = 500, number of meshes. High : +precise -fast / Low : -precise +fast
% windSpeed : (double) wind speed (km.h⁻¹)
% 
% example : EstimationRhopData_2part_Eul(400e-6, 900:1100, true, 600)

if nargin == 0
    SizePart = false; % particles size tested (m)
    RhoP_test = 800:100:1000; % plage de densite a tester (kg.m⁻³)
    saveFig = false;
    N = 500;
    windSpeed = 50;
end

%% Set Parameters
path = '../Results/EstimRho/'; % saved figures directory

% Multinet sample characteristics 
% Station = 'RN2';
% Date = datetime('3/18/2021');


%% Water column parameters
% ModeleHydro='2012RHOMA_arome_003.nc';
% SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
% % Load file with indexes corresponding to each stations
% stationFile = '../Data/stationIJ_CEREGE.mat';
% load(stationFile,'stationIJ');
% % Get the indexes of the right station
% I0 = stationIJ{stationIJ{:,'station'} == Station,'I0'};
% J0 = stationIJ{stationIJ{:,'station'} == Station,'J0'};
% load(SauvegardeModeleHydro, 'H0')
% L = H0(I0,J0); 
L = 50;
if nargin < 4
    N = 200; % number of meshes. High : +precise -fast / Low : -precise +fast
end
dz= L/N; % size of meshes
z = (0:dz:L)'; % meshes boundaries
z_ = z(1:end-1)+dz/2; % center of the meshes

% day = '3fev'; % day corresponding to diffusive turbulence data
day = '10fev';
% day = false;

% load('../Data/stationLonLat_CEREGE.mat','station') % load position data of stations
% staName = 'RN2'; % chosen station
% iStation = station(all(station(:,'Station').Variables'==staName')',:); % corresponding row of the table station
% % Position of station
% Lon0 = iStation(:,'Lon').Variables;
% Lat0 = iStation(:,'Lat').Variables;


%% Data

dh = 0.71; % Net oppening (m)

% SamplingDate = datetime('3/18/2021');
% DataFile = '../Data/data_mps.txt';
% [ConcentrationSample, DepthSample] = getDataNpart(false, SizePart, true, DataFile, SamplingDate);

% 3fev
% CMes=[0.27 0.08 0.09 0.1 0.2]; % Concentration
% ZMes=[0 25 35 45 50]+dh/2; % Depth of the samples

% 10fev
CMes=[0.62 0.34 0.06 0.02 0]; % Concentration
ZMes=[1 10 15 40 50]+dh/2; % Depth of the samples

ConcentrationSample = CMes(1:end)';
DepthSample = ZMes(1:end)';
% DataInterp = interp1(ZMes,CMes,z_,'pchip')'; % Interpolated data

% %% Data
% 
% dh = 0.71; % Net oppening (m)
% 
% % Concentration at each depth and filtered volume
% ConcVolFile = '../Data/ConcVol_MP.txt';
% ConcVolTable = load_ConcVol_data(ConcVolFile);
% 
% % Multinet condition
% concCond = strcmp(ConcVolTable{:,'station'},Station) & ConcVolTable{:,'date'} == Date;
% 
% % Concentration data
% CMes = ConcVolTable{concCond,'C'};
% ZMes = ConcVolTable{concCond,'depth'}+dh/2;
% 
% 
% ConcentrationSample = CMes(1:end);
% DepthSample = ZMes(1:end);
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


if SizePart
    modSize = SizePart;
else
    modSize = 350e-6;
end


%% Run model
clear Resultats

CRho = zeros(length(RhoP_test),length(ConcentrationSample)+1);
concRho = zeros(length(RhoP_test),N+1);
for i = 1:length(RhoP_test)
    rho=RhoP_test(i);
%     [conc, z_] = Transport_Eulerian(modSize, rho, N, L, day, windSpeed, Date);
    [conc, z_] = Transport_Eulerian(modSize, rho, N, L, day);
    
      % find the modeled concentration at depth corresponding with sample
    CModel = zeros(size(ConcentrationSample));
    j = 0;
    for ic = 1:2:length(ibound)
        j = j+1;
        CModel(j) = mean(conc(ibound(ic):ibound(ic+1)));
    end, clear ic j,
    
    CRho(i,:) = [rho CModel'];
    concRho(i,:) = [rho conc];
end, clear i,

Resultats(1:length(RhoP_test)) = struct('Rho1', 0, 'Rho2', 0,...
    'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'rmseErreur', 0);

combRho = [nchoosek(RhoP_test,2);[RhoP_test' RhoP_test']];

for iRes=1:size(combRho,1)
    Rho1 = combRho(iRes,1);
    Rho2 = combRho(iRes,2);
    
    Resultats(iRes).Rho1 = Rho1;
    Resultats(iRes).Rho2 = Rho2;
    Resultats(iRes).ConcentrationModel = (CRho(CRho(:,1)==Rho1,2:end) + CRho(CRho(:,1)==Rho2,2:end))';
    Resultats(iRes).conc = (concRho(concRho(:,1)==Rho1,2:end) + concRho(concRho(:,1)==Rho2,2:end))';
end, clear iRes,

%% Model output treatment
for i=1:length(Resultats)

    alpha = Resultats(i).ConcentrationModel\ConcentrationSample;
%     alpha = Resultats(i).conc\DataInterp';

    Erreur = abs(Resultats(i).ConcentrationModel.*alpha - ConcentrationSample);
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
    
    Resultats(i).Alpha = alpha;
    Resultats(i).Erreur = Erreur;
    Resultats(i).rmseErreur = rmseErreur;

end, clear i,

ErrorPlotRes = zeros(length(RhoP_test));
for i=1:length(Resultats)
    a = RhoP_test == Resultats(i).Rho1;
    b = RhoP_test == Resultats(i).Rho2;
    ErrorPlotRes(a,b) = Resultats(i).rmseErreur;
end, clear i,

ErrorPlotRes = ErrorPlotRes + ErrorPlotRes' - eye(size(ErrorPlotRes)).*ErrorPlotRes;

ErrorPlotResBis = ErrorPlotRes;
% if nargin < 6
%     ErrorPlotResBis = ErrorPlotRes;
% else
%     ErrorPlotResBis = NaN(length(RhoP_test));
%     ErrorPlotResBis(1:(length(sub1)-1),(end-length(sub2)+1):end) = ErrorPlotRes(1:(length(sub1)-1),(end-length(sub2)+1):end);
%     ErrorPlotResBis((end-length(sub2)+1):end,1:(length(sub1)-1)) = ErrorPlotRes((end-length(sub2)+1):end,1:(length(sub1)-1));
% end


minI = [Resultats.rmseErreur] == min([Resultats.rmseErreur]);
minRho1 = Resultats(minI).Rho1;
minRho2 = Resultats(minI).Rho2;


%% Display results
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm'];

% f1 = figure(1); clf,
% plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
% xlim([0 4])
% ylim([-L 0])
% xlabel('Concentration (mps.m⁻¹)')
% ylabel('Depth (m)')
% hold on
% plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
% for res = Resultats
%     plot(res.conc*res.Alpha,-z_,'DisplayName',  ['Rho1 = ' num2str(res.Rho1) ', Rho2 = ' num2str(res.Rho2) 'kg.m⁻³'])
% end
% legend('Location', 'best')
% title(ttl)
% hold off

% figure 2 : plot error with respect to rhoP
f2 = figure(2); clf,
pcolor(RhoP_test,RhoP_test,ErrorPlotResBis);
c = colorbar;
c.Label.String = 'RMSE (mps.m⁻³)';
xlabel('Rho_1 (kg.m⁻³)');
ylabel('Rho_2 (kg.m⁻³)');
% xticks(RhoP_test);
% yticks(RhoP_test);
title(ttl)

% caxis([0.01 0.13])

% figure 3 : plot best rhoP profile
f3 = figure(3); clf,
% plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
xlim([0 7])
ylim([-L 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on 
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
%plot(Resultats(minI).ConcentrationModel*Resultats(minI).Alpha,-DepthSample,'DisplayName', ['Rho1 = ' num2str(minRho1) ', Rho2 = ' num2str(minRho2) 'kg.m⁻³'])
plot(Resultats(minI).conc*Resultats(minI).Alpha,-z_,'DisplayName',  ['Rho1 = ' num2str(minRho1) ', Rho2 = ' num2str(minRho2) 'kg.m⁻³'])
legend('Location', 'best')
title([ttl ' -- Rho1 = ' num2str(minRho1) ', Rho2 = ' num2str(minRho2) 'kg.m⁻³'])
hold off

%% Save figures
if saveFig
    if length(RhoP_test) == 1
        rhoInter = num2str(RhoP_test);
    else
        rhoInter = [num2str(min(RhoP_test)) '-' num2str(RhoP_test(2)-RhoP_test(1)) '-' num2str(max(RhoP_test))];
    end
    
    fileName = ['size' num2str(modSize*1e6) '_rho' rhoInter '_data10FevMarine'];
    
    F = [f2, f3];
    xPart = '2part_';
    N = {'error_', 'min_'};
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
