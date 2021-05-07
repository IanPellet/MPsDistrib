function [Resultats,ConcentrationSample,Ztest_,z_] = EstimationRhopData_2part_Eul(SizePart, type_name, RhoP_test,saveFig, sub1, sub2)

%% Set Parameters

path = '../Results/EstimRho/';

% Wind parameters
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

load('../Data/stationLonLat_CEREGE.mat','station') % load position data of stations
staName = 'RN2'; % chosen station
iStation = station(all(station(:,'Station').Variables'==staName')',:); % corresponding row of the table station
% Position of station
Lon0 = iStation(:,'Lon').Variables;
Lat0 = iStation(:,'Lat').Variables;


%% Data loading

L = 50;
N = 500;
dz= L/N;
z = 0:dz:L;
z_ = z(1:end-1)+dz/2;
dh = 0.15; % Net oppening (m)

% [CMes, ZMes] = getDataNpart(type_name, SizePart, true);
% CMes = CMes(1:end-1);
% ZMes = ZMes(1:end-1);

CMes=[0.27 0.08 0.09 0.1 0.2];
ZMes=[0 25 35 45 50]+dh;
ConcentrationSample = interp1(ZMes,CMes,z_,'pchip')';
DepthSample = z_';
day = '3fev';



% Boundaries in between which modeled part number have to be tested
boundTest = zeros(1,length(DepthSample)*2);
for i = 1:length(boundTest)
    if mod(i,2)
        boundTest(i) = DepthSample(fix((i+1)/2))-dh;
    else
        boundTest(i) = DepthSample(fix((i+1)/2));
    end
end, clear i,

Ztest_ = DepthSample - dh/2;


if SizePart
    modSize = SizePart;
else
    modSize = 350e-6;
end


%% Run model
clear Resultats

CRho = zeros(length(RhoP_test),N+1);
for i = 1:length(RhoP_test)
    rho=RhoP_test(i);
    [C, z_] = Transport_Eulerian(modSize, rho, N, day);
    CRho(i,:) = [rho C];
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
end, clear iRes,

%% Model output treatment
for i=1:length(Resultats)

    alpha = Resultats(i).ConcentrationModel\ConcentrationSample;
    
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

% for i=1:length(RhoP_test)
%     for j=1:length(RhoP_test)
%         r1 = RhoP_test(i);
%         r2 = RhoP_test(j);
%         
%         iCond = find(([Resultats.Rho1] == r1 & [Resultats.Rho2] == r2) |...
%             ([Resultats.Rho1] == r2 & [Resultats.Rho2] == r1),1);
%         
%         ErrorPlotRes(i,j) = Resultats(iCond).rmseErreur;
%     end, clear j,
% end, clear i,
%               

%% Display results
% if type_name
%     dispType = type_name;
% else
%     dispType = 'all';
% end
% ttl = ['Particle size : ' num2str(modSize*1e6) 'µm -- Particle type : ' dispType];
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm'];

f1 = figure(1); clf,
plot(ConcentrationSample,-Ztest_,'--', 'DisplayName', 'Data interpolation');
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
for res = Resultats
    plot(res.ConcentrationModel*res.Alpha,-z_,'DisplayName', ['Rho1 = ' num2str(res.Rho1) ', Rho2 = ' num2str(res.Rho2) 'kg.m⁻³'])
end
legend('Location', 'best')
title(ttl)
hold off

if nargin < 6
    ErrorPlotResBis = ErrorPlotRes;
else
    ErrorPlotResBis = NaN(length(RhoP_test));
    ErrorPlotResBis(1:(length(sub1)-1),(end-length(sub2)+1):end) = ErrorPlotRes(1:(length(sub1)-1),(end-length(sub2)+1):end);
    ErrorPlotResBis((end-length(sub2)+1):end,1:(length(sub1)-1)) = ErrorPlotRes((end-length(sub2)+1):end,1:(length(sub1)-1));
end

f2 = figure(2); clf,
pcolor(RhoP_test,RhoP_test,ErrorPlotResBis);
c = colorbar;
c.Label.String = 'RMSE (mps.m⁻³)';
xlabel('Rho_1 (kg.m⁻³)');
ylabel('Rho_2 (kg.m⁻³)');
xticks(RhoP_test);
yticks(RhoP_test);


title(ttl)


minI = [Resultats.rmseErreur] == min([Resultats.rmseErreur]);
minRho1 = Resultats(minI).Rho1;
minRho2 = Resultats(minI).Rho2;


f3 = figure(3); clf,
plot(ConcentrationSample,-Ztest_,'--', 'DisplayName', 'Data interpolation');
% xlim([0 max(ConcentrationSample)])
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on 
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
plot(Resultats(minI).ConcentrationModel*Resultats(minI).Alpha,-z_,'DisplayName', ['Rho1 = ' num2str(minRho1) ', Rho2 = ' num2str(minRho2) 'kg.m⁻³'])
legend('Location', 'best')
title([ttl ' -- Rho1 = ' num2str(minRho1) ', Rho2 = ' num2str(minRho2) 'kg.m⁻³'])
hold off

if saveFig
    if length(RhoP_test) == 1
        rhoInter = num2str(RhoP_test);
    else
        rhoInter = [num2str(min(RhoP_test)) '-' num2str(RhoP_test(2)-RhoP_test(1)) '-' num2str(max(RhoP_test))];
    end
    
    fileName = ['size' num2str(modSize*1e6) '_rho' rhoInter 'part_EulMarine'];
    
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
