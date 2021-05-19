function [Resultats,ConcentrationSample,Ztest_,z_] = EstimationRhopData_1part_Lagr(SizePart, type_name, RhoP_test, saveFig, N)

%% Set Parameters

path = '../Results/EstimRho/';

wind = 1; % km/h
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

nPart = 50e3; % number of particles

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
% N = 100;
dz= L/N;
z = (0:dz:L)'; % meshes boundaries
z_ = z(1:end-1)+dz/2; % center of the meshes

day = '10fev'; % day corresponding to diffusive turbulence data

% [ConcentrationSample, DepthSample] = getDataNpart(type_name, SizePart, true);
CMes=[0.62 0.34 0.06 0.02 0]; % Concentration
ZMes=[1 10 15 40 50]; % Depth of the samples
ConcentrationSample = CMes';
DepthSample = ZMes';
DataInterp = interp1(ZMes,CMes,z_,'pchip')'; % Interpolated data

dh = 0.71; % Net oppening (m)

%% Find corresponding depths to get from the model
boundTest = zeros(length(DepthSample)*2,1); % Boundaries in between which modeled part number have to be tested
ibound = zeros(size(boundTest)); % index boundaries of meshes to avg

for i = 1:length(boundTest)
    if mod(i,2)
        boundTest(i) = DepthSample(fix((i+1)/2));
    else
        boundTest(i) = DepthSample(fix((i+1)/2))+dh;
    end
    ibound(i) = min(length(z_),max(1,fix(boundTest(i)/dz)+1));
end, clear i,

Ztest_ = DepthSample + dh/2; 

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

NsecTest = 60*60;
% launch simulation for each rhop
for iRes = 1:length(RhoP_test)
    RhoP = RhoP_test(iRes);    
    [~, ~, ~, ppHist] = varMP_model(modSize, RhoP, TypePart, nPart, tf, dt_test, wind, month, Lon0, Lat0,L, N, day,NsecTest);
    
    hConc = NaN(size(ppHist,1),length(z_));
    for hStep = 1:size(ppHist,1)
        pp = ppHist(hStep,:);
        if ~isnan(sum(pp))
            histi = histogram(pp, "BinEdges", z, 'Visible', 'off').Values;
            hConc(hStep,:) = histi/dz*L/nPart;
        end
    end, clear hStep,

    conc = mean(hConc, 'omitnan');
    stdConc = std(hConc, 'omitnan');
    
    ConcentrationModel = zeros(size(ConcentrationSample));
    j = 0;
    for ic = 1:2:length(ibound)
        j = j+1;
        ConcentrationModel(j) = mean(conc(ibound(ic):ibound(ic+1)));
    end, clear ic j,
    
    
%     % 2) calculer l'erreur avec les mesures :
%     % =======================================
%     % find the modeled concentration at depth corresponding with sample
%     hTest = histogram(PartPos, "BinEdges",  boundTest, 'Visible', 'off').Values;
%     NpartModel = zeros(size(ConcentrationSample));
%     j = 0;
%     for i = 1:length(hTest)
%         if mod(i,2)
%             j = j+1;
%             NpartModel(j) = hTest(i);
%         end
%     end
% %     conc = histogram(PartPos, "BinEdges",z,'Visible', 'off').Values/dz / (nPart/L) ;
%     
%     ConcentrationModel = NpartModel/dh / (nPart/L) ; % on ramène le modèle à 1
% %     alpho = mean(ConcentrationSample(ConcentrationModel>0)./ConcentrationModel(ConcentrationModel>0));
% %     alpha = ConcentrationSample\ConcentrationModel;
    alpha = ConcentrationModel\ConcentrationSample;
    
    
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
%     rmseErreur = mean(Erreur,'omitnan');
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
    
    Resultats(iRes).RhoP = RhoP;
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).Erreur = Erreur;
    Resultats(iRes).rmseErreur = rmseErreur;
    Resultats(iRes).conc = conc;
    Resultats(iRes).std = stdConc;
end, clear iRes RhoP,

% Find rhop correponding to minimal error
minI = [Resultats.rmseErreur] == min([Resultats.rmseErreur]); 
minRho = Resultats(minI).RhoP;

%% Display results
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm']; % Figures titles

% figure 1 : plot all simulation results
f1 = figure(1); clf,
plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
xlim([0 max(Resultats(minI).conc*Resultats(minI).Alpha)])
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
for res = Resultats
    plot(res.conc*res.Alpha,-z_,'DisplayName', ['RhoP = ' num2str(res.RhoP) 'kg.m⁻³'])
end
plot((Resultats(minI).conc+2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', 'min+2std');
plot((Resultats(minI).conc-2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', 'min-2std');
legend('Location', 'southeast')
title(ttl)
hold off

% f1 = figure(1); clf,
% plot((Resultats(minI).conc+Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--r',(Resultats(minI).conc-Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--r');
% xlim([0 max(Resultats(minI).conc*Resultats(minI).Alpha)])
% ylim([-L+0.75 0])
% xlabel('Concentration (mps.m⁻¹)')
% ylabel('Depth (m)')
% hold on
% for res = Resultats
%     plot(res.conc*res.Alpha,-z_,'DisplayName', ['RhoP = ' num2str(res.RhoP) 'kg.m⁻³'])
% end
% legend('Location', 'southeast')
% title(ttl)
% hold off

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
plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
xlim([0 max(Resultats(minI).conc*Resultats(minI).Alpha)])
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on 
plot(CMes,-Ztest_,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
plot(Resultats(minI).conc*Resultats(minI).Alpha,-z_,'DisplayName', ['RhoP = ' num2str(Resultats(minI).RhoP) 'kg.m⁻³'])
plot((Resultats(minI).conc+2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', '+2std');
plot((Resultats(minI).conc-2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', '-2std');
legend('Location', 'southeast')
title([ttl ' -- Rho_p = ' num2str(minRho) 'kg.m⁻³'])
hold off

%% Save figures
if saveFig
    if length(RhoP_test) == 1
        rhoInter = num2str(RhoP_test);
    else
        rhoInter = [num2str(min(RhoP_test)) '-' num2str(RhoP_test(2)-RhoP_test(1)) '-' num2str(max(RhoP_test))];
    end
    
    fileName = ['size' num2str(modSize*1e6) '_rho' rhoInter '_lagr'];
    
    F = [f1, f2, f3];
    xPart = '1part_';
    N = {'profils_', 'error_', 'min_'};Dvar
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
