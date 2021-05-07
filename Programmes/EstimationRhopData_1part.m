function [Resultats,ConcentrationSample,Ztest_,z_] = EstimationRhopData_1part(SizePart, type_name, RhoP_test, saveFig)

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
z = 0:dz:L;
z_ = z(1:end-1)+dz/2;


% [ConcentrationSample, DepthSample] = getDataNpart(type_name, SizePart, true);
CMes=[0.62 0.34 0.06 0.02 0];
ZMes=[1 10 15 40 L];
ConcentrationSample = interp1(ZMes,CMes,z_,'pchip')';
DepthSample = z_';

dh = 0.15; % Net oppening (m)

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

clear Resultats

Resultats(1:length(RhoP_test)) = struct('RhoP', 0, 'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'rmseErreur', 0);
if SizePart
    modSize = SizePart;
else
    modSize = 350e-6;
end


iRes = 0;
for RhoP = RhoP_test
    iRes = iRes+1;
    [~, ~, PartPos] = varMP_model(modSize, RhoP, TypePart, nPart, tf, dt_test, wind, month, Lon0, Lat0,L, N, path);
%     
%     figure(2)
%     hist = histogram(PartPos,'Visible', 'off').Values;
%     plot(hist,-linspace(0,L,length(hist)))
    
    % 2) calculer l'erreur avec les mesures :
    % =======================================
    hTest = histogram(PartPos, "BinEdges",  boundTest, 'Visible', 'off').Values;
%     NpartModel = hcalc(fix(DepthSample./dh));
    NpartModel = zeros(size(ConcentrationSample));
    j = 0;
    for i = 1:length(hTest)
        if mod(i,2)
            j = j+1;
            NpartModel(j) = hTest(i);
        end
    end
    conc = histogram(PartPos, "BinEdges",z,'Visible', 'off').Values/dz * (nPart/L) ;
    
    ConcentrationModel = NpartModel/dh * (nPart/L) ; % on ramène le modèle à 1
%     alpho = mean(ConcentrationSample(ConcentrationModel>0)./ConcentrationModel(ConcentrationModel>0));
%     alpha = ConcentrationSample\ConcentrationModel;
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
end

% typeSTR = {'fibre','fragment','film','mousse'}; % List of existing type
% if type_name
%     dispType = type_name;
% else
%     dispType = 'all';
% end
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm'];

f1 = figure(1); clf,
plot(ConcentrationSample,-Ztest_,'--', 'DisplayName', 'Data interpolation');
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
for res = Resultats
    plot(res.conc*res.Alpha,-z_,'DisplayName', ['RhoP = ' num2str(res.RhoP) 'kg.m⁻³'])
%     plot(res.ConcentrationModel*res.Alpha, -Ztest_,'*','DisplayName', ['RhoP = ' num2str(res.RhoP)])
end
legend('Location', 'southeast')
title(ttl)
hold off

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

minI = [Resultats.rmseErreur] == min([Resultats.rmseErreur]);
minRho = Resultats(minI).RhoP;


f3 = figure(3); clf,
plot(ConcentrationSample,-Ztest_,'--', 'DisplayName', 'Data interpolation');
% xlim([0 max(ConcentrationSample)])
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on 
plot(CMes,-ZMes,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
plot(Resultats(minI).conc*Resultats(minI).Alpha,-z_,'DisplayName', ['RhoP = ' num2str(Resultats(minI).RhoP) 'kg.m⁻³'])
legend('Location', 'southeast')
title([ttl ' -- Rho_p = ' num2str(minRho) 'kg.m⁻³'])
hold off

if saveFig
    if length(RhoP_test) == 1
        rhoInter = num2str(RhoP_test);
    else
        rhoInter = [num2str(min(RhoP_test)) '-' num2str(RhoP_test(2)-RhoP_test(1)) '-' num2str(max(RhoP_test))];
    end
    
    fileName = ['size' num2str(modSize*1e6) '_rho' rhoInter 'part_10fev'];
    
    F = [f1, f2, f3];
    xPart = '1part_';
    N = {'profils_', 'error_', 'min_'};
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
