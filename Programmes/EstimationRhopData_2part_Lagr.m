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
L = 55;
N = L;
dz= L/N;
z = (0:dz:L)'; % meshes boundaries
z_ = z(1:end-1)+dz/2; % center of the meshes

day = '3fev'; % day corresponding to diffusive turbulence data

%% Data
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
NsecTest = 60*30;
ppRho(length(RhoP_test)) = struct('rho', 0);
for i = 1:length(RhoP_test)
    rho=RhoP_test(i);
    [~, ~, ~, ppHist] = varMP_model(modSize, rho, TypePart, nPart/2, tf, dt_test, wind, month, Lon0, Lat0,L, N, day,NsecTest);
    ppRho(i).rho = rho;
    ppRho(i).ppHist = ppHist;
end, clear i,

Resultats(1:length(RhoP_test)) = struct('Rho1', 0, 'Rho2', 0,...
    'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'rmseErreur', 0);


%% Model output treatment
combRho = [nchoosek(RhoP_test,2);[RhoP_test' RhoP_test']];

for iRes=1:size(combRho,1)
    Rho1 = combRho(iRes,1);
    Rho2 = combRho(iRes,2);    
    Resultats(iRes).Rho1 = Rho1;
    Resultats(iRes).Rho2 = Rho2;
    
    pphist = [ppRho([ppRho(:).rho] == Rho1).ppHist  ppRho([ppRho(:).rho] == Rho2).ppHist];
    nPP = size(pphist,1);
    pphist = reshape(pphist, [numel(pphist) 1]);
%     hConc = NaN(size(pphist,1),length(z_));
%     for hStep = 1:size(ppHist,1)
%         pp = pphist(hStep,:);
%         [histi,~] =  groupcounts(pp',z,'IncludeEmptyGroups',true);
%         hConc(hStep,:) = histi'/dz*L/nPart;
%     end, clear hStep pp histi,
    
    [histi,~] =  groupcounts(pphist,z,'IncludeEmptyGroups',true);
    hConc = histi'/dz*L/nPart;
    conc = hConc/nPP;
%     conc = mean(hConc, 'omitnan');
%     stdConc = std(hConc, 'omitnan');

    ConcentrationModel = zeros(size(ConcentrationSample));
    j = 0;
    for ic = 1:2:length(ibound)
        j = j+1;
        ConcentrationModel(j) = mean(conc(ibound(ic):ibound(ic+1)));
    end, clear ic j,

    alpha = ConcentrationModel\ConcentrationSample;
    % Compute error
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));

    % Fill resultats
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).rmseErreur = rmseErreur;
    Resultats(iRes).conc = conc;
%     Resultats(iRes).std = stdConc;

end, clear iRes pphist hConc hConc conc stdConc ConcentrationModel nPP,

%% Create error matrix to plot results

ErrorPlotRes = zeros(length(RhoP_test));
for i=1:length(Resultats)
    a = RhoP_test == Resultats(i).Rho1;
    b = RhoP_test == Resultats(i).Rho2;
    ErrorPlotRes(a,b) = Resultats(i).rmseErreur;
end, clear i a b,

ErrorPlotResBis = ErrorPlotRes + ErrorPlotRes' - eye(size(ErrorPlotRes)).*ErrorPlotRes;
        


minI = find([Resultats.rmseErreur] == min([Resultats.rmseErreur]),1);

minRho1 = Resultats(minI).Rho1;
minRho2 = Resultats(minI).Rho2;

%% Display results
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm'];

% f1 = figure(1); clf,
% plot(DataInterp,-z_,'--', 'DisplayName', 'Data interpolation');
% ylim([-L+0.75 0])
% xlabel('Concentration (mps.m⁻¹)')
% ylabel('Depth (m)')
% hold on
% plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
% for res = Resultats
%     plot(res.conc*res.Alpha,-z_,'DisplayName', ['Rho1 = ' num2str(res.Rho1) ', Rho2 = ' num2str(res.Rho2) 'kg.m⁻³'])
% end
% plot((Resultats(minI).conc+2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', 'min+2std');
% plot((Resultats(minI).conc-2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', 'min-2std');
% legend('Location', 'southeast')
% title(ttl)
% hold off

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
% plot((Resultats(minI).conc+2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', '+2std');
% plot((Resultats(minI).conc-2*Resultats(minI).std)*Resultats(minI).Alpha,-z_, '--','DisplayName', '-2std');
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
