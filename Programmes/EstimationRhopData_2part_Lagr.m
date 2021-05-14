function [Resultats,ConcentrationSample,Ztest_,z_] = EstimationRhopData_2part_Lagr(SizePart, type_name, RhoP_test,saveFig)

%% Set Parameters

path = '../Results/EstimRho/';

wind = 50; % km/h
month = 03; 


if nargin == 0
    SizePart = false; % particles size tested (m)
    type_name = false;
    RhoP_test = 800:100:1000; % plage de densite a tester (kg.m⁻³)
    saveFig = false;
end

if type_name
    type_dict = containers.Map({'fibre','fragment','film','mousse'},0:3); % possible types
    TypePart = type_dict(type_name); % corresponding int
else
    TypePart = false;
end


nPart = 10e3; % number of particles

tf = 1e5; % simulation time (s)
dt_test = 60*60; % test time interval (s)

load('../Data/stationLonLat_CEREGE.mat','station') % load position data of stations
staName = 'RN2'; % chosen station
iStation = station(all(station(:,'Station').Variables'==staName')',:); % corresponding row of the table station
% Position of station
Lon0 = iStation(:,'Lon').Variables;
Lat0 = iStation(:,'Lat').Variables;


%% Data loading

% ModeleHydro='2012RHOMA_arome_003.nc';
% SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
% load(SauvegardeModeleHydro)
% [I0,J0]=ReperePoint(Lon,Lat,Lon0,Lat0);
% L = H0(I0,J0);
L = 55;
N = L;
dz= L/N;
z = (0:dz:L)'; % meshes boundaries
z_ = z(1:end-1)+dz/2; % center of the meshes

day = '3fev'; % day corresponding to diffusive turbulence data

%% Data

% SamplingDate = datetime('3/18/2021');
% DataFile = '../Data/data_mps.txt';
% [ConcentrationSample, DepthSample] = getDataNpart(false, SizePart, true, DataFile, SamplingDate);
CMes=[0.27 0.08 0.09 0.1 0.2]; % Concentration
ZMes=[0 25 35 45 50]; % Depth of the samples
ConcentrationSample = CMes(1:end)';
DepthSample = ZMes(1:end)';
DataInterp = interp1(ZMes,CMes,z_,'pchip')'; % Interpolated data

dh = 0.71; % Net oppening (m)


if SizePart
    modSize = SizePart;
else
    modSize = 350e-6;
end


%% Find corresponding depths to get from the model
boundTest = zeros(length(DepthSample)*2,1); % Boundaries in between which modeled part number have to be tested
% ibound = zeros(size(boundTest)); % index boundaries of meshes to avg

for i = 1:length(boundTest)
    if mod(i,2)
        boundTest(i) = DepthSample(fix((i+1)/2));
    else
        boundTest(i) = DepthSample(fix((i+1)/2))+dh;
    end
%     ibound(i) = min(length(z_),max(1,fix(boundTest(i)/dz)+1));
end, clear i,

Ztest_ = DepthSample + dh/2; 

%% Simulations
clear Resultats

pPosRho = zeros(length(RhoP_test),nPart/2+1);
for i = 1:length(RhoP_test)
    rho=RhoP_test(i);
    [~, ~, pPos] = varMP_model(modSize, rho, TypePart, nPart/2, tf, dt_test, wind, month, Lon0, Lat0, L, N, day);
    pPosRho(i,:) = [rho pPos];
end, clear i,

Resultats(1:length(RhoP_test)) = struct('Rho1', 0, 'Rho2', 0, 'PartPos',...
    [] ,'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'rmseErreur', 0);


%% Model output treatment
combRho = [nchoosek(RhoP_test,2);[RhoP_test' RhoP_test']];

for iRes=1:size(combRho,1)
    
    Rho1 = combRho(iRes,1);
    Rho2 = combRho(iRes,2);    
    Resultats(iRes).Rho1 = Rho1;
    Resultats(iRes).Rho2 = Rho2;
    Resultats(iRes).PartPos = [pPosRho(pPosRho(:,1)==Rho1,2:end) pPosRho(pPosRho(:,1)==Rho2,2:end)];
    
    % Compute concentrations   
    hTest = histogram(Resultats(iRes).PartPos, "BinEdges", boundTest,'Visible', 'off').Values';
    NpartModel = zeros(size(ConcentrationSample));
    j = 0;
    for i = 1:length(hTest)
        if mod(i,2)
            j = j+1;
            NpartModel(j) = hTest(i);
        end
    end
    ConcentrationModel = NpartModel/dh;
    conc = histogram(Resultats(iRes).PartPos, "BinEdges",z,'Visible', 'off').Values/dz * (nPart/L) ;
    
    % Find best proportionnality coef between model output and data
%     Cmod = ConcentrationModel;
%     Cmod(Cmod == 0) = 1e-100;
%     alpha = ConcentrationModel\ConcentrationSample;
    alpha = conc'\DataInterp';
    % Compute error
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
    
    % Fill resultates
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).Erreur = Erreur;
    Resultats(iRes).rmseErreur = rmseErreur;
    Resultats(iRes).conc = conc;
end, clear iRes,

%% Create error matrix to plot results

ErrorPlotRes = zeros(length(RhoP_test));
for i=1:length(Resultats)
    a = RhoP_test == Resultats(i).Rho1;
    b = RhoP_test == Resultats(i).Rho2;
    ErrorPlotRes(a,b) = Resultats(i).rmseErreur;
end, clear i,

ErrorPlotResBis = ErrorPlotRes + ErrorPlotRes' - eye(size(ErrorPlotRes)).*ErrorPlotRes;
        

%% Display results
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm'];

f1 = figure(1); clf,
plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
for res = Resultats
    plot(res.conc*res.Alpha,-z_,'DisplayName', ['Rho1 = ' num2str(res.Rho1) ', Rho2 = ' num2str(res.Rho2) 'kg.m⁻³'])
end
legend('Location', 'southeast')
title(ttl)
hold off

% figure 2 : plot error with respect to rhoP
f2 = figure(2); clf,
pcolor(RhoP_test,RhoP_test,ErrorPlotResBis);
c = colorbar;
c.Label.String = 'RMSE (mps.m⁻³)';
xlabel('Rho_1 (kg.m⁻³)');
ylabel('Rho_2 (kg.m⁻³)');
xticks(RhoP_test);
yticks(RhoP_test);
title(ttl)


minI = find([Resultats.rmseErreur] == min([Resultats.rmseErreur]),1);

minRho1 = Resultats(minI).Rho1;
minRho2 = Resultats(minI).Rho2;

% figure 3 : plot best rhoP profile
f3 = figure(3); clf,
plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
% xlim([0 max(ConcentrationSample)])
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on 
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
plot(Resultats(minI).conc*Resultats(minI).Alpha,-z_,'DisplayName', ['Rho1 = ' num2str(minRho1) ', Rho2 = ' num2str(minRho2) 'kg.m⁻³'])
legend('Location', 'southeast')
title([ttl ' -- Rho1 = ' num2str(minRho1) ', Rho2 = ' num2str(minRho2) 'kg.m⁻³'])
hold off

if saveFig
    if length(RhoP_test) == 1
        rhoInter = num2str(RhoP_test);
    else
        rhoInter = [num2str(min(RhoP_test)) '-' num2str(RhoP_test(2)-RhoP_test(1)) '-' num2str(max(RhoP_test))];
    end
    
    fileName = ['size' num2str(modSize*1e6) '_rho' rhoInter '_MarineLagr'];
    
    F = [f2, f3];
    xPart = '2part_';
    N = {'error_', 'min_'};
    for i=1:length(F)
        f = F(i);
        n = N{i};
        NamePNG = [path 'png/' xPart n fileName '.png'];
        NameEPS = [path 'eps/' xPart n fileName '.eps'];
        NameFIG = [path 'fig/' xPart n fileName '.fig'];
        exportgraphics(f,NamePNG);
        exportgraphics(f,NameEPS,'ContentType','vector');
        savefig(f,NameFIG);
    end, clear i,
end

end
