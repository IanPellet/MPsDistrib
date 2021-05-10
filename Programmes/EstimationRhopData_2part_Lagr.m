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
L = 51;
N = L;
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
% for i = 1:length(boundTest)
%     temp_depth = DepthSample(fix((i+1)/2));
%     if i>1 && boundTest(i-1) == temp_depth
%         boundTest(i) = temp_depth+dh;
%     else
%         boundTest(i) = temp_depth;
%     end
% end
% boundTest = 0:dh:L;
for i = 1:length(boundTest)
    if mod(i,2)
        boundTest(i) = DepthSample(fix((i+1)/2))-dh;
    else
        boundTest(i) = DepthSample(fix((i+1)/2));
    end
end, clear i,

% Ztest_ = DepthSample + dh/2;
Ztest_ = DepthSample - dh/2;

nPart = 50e3; % number of particles

tf = 1e5; % simulation time (s)
dt_test = 60*60; % test time interval (s)


if SizePart
    modSize = SizePart;
else
    modSize = 350e-6;
end


%% Run model
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
% 
% combRho = nchoosek(RhoP_test,2);



%% Model output treatment
combRho = [nchoosek(RhoP_test,2);[RhoP_test' RhoP_test']];

for iRes=1:size(combRho,1)
    
    Rho1 = combRho(iRes,1);
    Rho2 = combRho(iRes,2);    
    Resultats(iRes).Rho1 = Rho1;
    Resultats(iRes).Rho2 = Rho2;
    Resultats(iRes).PartPos = [pPosRho(pPosRho(:,1)==Rho1,2:end) pPosRho(pPosRho(:,1)==Rho2,2:end)];
    
    % Compute concentrations   
    ConcentrationModel = histogram(Resultats(iRes).PartPos, "BinEdges",z,'Visible', 'off').Values'/dz;
    % Find best proportionnality coef between model output and data
%     Cmod = ConcentrationModel;
%     Cmod(Cmod == 0) = 1e-100;
    alpha = ConcentrationModel\ConcentrationSample;
    % Compute error
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
    
    % Fill resultates
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).Erreur = Erreur;
    Resultats(iRes).rmseErreur = rmseErreur;
    
end, clear iRes,

%% Create error matrix to plot results

ErrorPlotRes = zeros(length(RhoP_test));
for i=1:length(Resultats)
    a = RhoP_test == Resultats(i).Rho1;
    b = RhoP_test == Resultats(i).Rho2;
    ErrorPlotRes(a,b) = Resultats(i).rmseErreur;
end, clear i,

ErrorPlotResBis = ErrorPlotRes + ErrorPlotRes' - eye(size(ErrorPlotRes)).*ErrorPlotRes;

% ErrorPlotRes = zeros(length(RhoP_test));
% for i=1:length(RhoP_test)
%     for j=1:length(RhoP_test)
%         r1 = RhoP_test(i);
%         r2 = RhoP_test(j);
%         partpos = [pPosRho(pPosRho(:,1)==r1,2:end) pPosRho(pPosRho(:,1)==r2,2:end)];
%         
%         hTest = histogram(partpos, "BinEdges",  boundTest, 'Visible', 'off').Values;
%         NpartModel = zeros(size(ConcentrationSample));
%         k = 0;
%         for in = 1:length(hTest)
%             if mod(in,2)
%                 k = k+1;
%                 NpartModel(k) = hTest(in);
%             end
%         end, clear k in,
% 
%         Cmod = NpartModel/dh * (nPart/L) ; % on ramène le modèle à 1
%         alpha = Cmod\ConcentrationSample;
% 
%         erreur = abs(Cmod.*alpha - ConcentrationSample);
%         rmse = sqrt(mean(erreur.^2,'omitnan'));
%         
%         ErrorPlotRes(i,j) = rmse;
%     end, clear j,
% end, clear i,
%         
% 
% for i=1:length(Resultats)
%     hTest = histogram(Resultats(i).PartPos, "BinEdges",  boundTest, 'Visible', 'off').Values;
%     NpartModel = zeros(size(ConcentrationSample));
%     j = 0;
%     for in = 1:length(hTest)
%         if mod(in,2)
%             j = j+1;
%             NpartModel(j) = hTest(in);
%         end
%     end
%     conc = histogram(Resultats(i).PartPos, "BinEdges",z,'Visible', 'off').Values/dz * (nPart/L) ;
%     
%     ConcentrationModel = NpartModel/dh * (nPart/L) ; % on ramène le modèle à 1
%     alpha = ConcentrationModel\ConcentrationSample;
%     
%     
%     Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
% %     rmseErreur = mean(Erreur,'omitnan');
%     rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
%     
%     Resultats(i).ConcentrationModel = ConcentrationModel;
%     Resultats(i).Alpha = alpha;
%     Resultats(i).Erreur = Erreur;
%     Resultats(i).rmseErreur = rmseErreur;
%     Resultats(i).conc = ConcentrationModel;
% end, clear i,


        

%% Display results
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm'];

% f1 = figure(1); clf,
% plot(ConcentrationSample,-Ztest_,'--', 'DisplayName', 'Data interpolation');
% ylim([-L+0.75 0])
% xlabel('Concentration (mps.m⁻¹)')
% ylabel('Depth (m)')
% hold on
% plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
% for res = Resultats
%     plot(res.ConcentrationModel*res.Alpha,-z_,'DisplayName', ['Rho1 = ' num2str(res.Rho1) ', Rho2 = ' num2str(res.Rho2) 'kg.m⁻³'])
% end
% legend('Location', 'southeast')
% title(ttl)
% hold off

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
