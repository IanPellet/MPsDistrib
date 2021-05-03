path = '../Results/EstimRho/';

wind = 1; % km/h
month = 03; 

load('../Data/stationLonLat_CEREGE.mat','station') % load position data of stations
staName = 'RN2'; % chosen station
iStation = find(station(:,'Station').Variables == staName); % corresponding row of the table station
% Position of station
Lon0 = station(1,'Lon').Variables;
Lat0 = station(1,'Lat').Variables;

ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)


% [I0,J0]=ReperePoint(Lon,Lat,Lon0,Lat0);
% L = H0(I0,J0);
L = 55;

[ConcentrationSample, DepthSample] = getDataConcentration;
% DepthSample = [0 5 15 25 35 55]; % depth of samples (m)
% ConcentrationSample = [0.99117913922709/0.75 ,1.16012084592145,... % le manta à 0m ne récupère que 75% des particules
%     0.937931034482759,0.567505720823799,1.3265306122449,1.18213058419244]; % mps.m⁻³
% FilteredVolume = [88.78314375,331,435,437,392,291]; % m³


dh = 0.15; % Net oppening (m)

% Boundaries in between which modeled part number have to be tested
boundTest = zeros(1,length(DepthSample)*2);
for i = 1:length(boundTest)
    temp_depth = DepthSample(fix((i+1)/2));
    if i>1 && boundTest(i-1) == temp_depth
        boundTest(i) = temp_depth+dh;
    else
        boundTest(i) = temp_depth;
    end
end

Ztest_ = DepthSample + dh/2;

SizePart = 750*1e-6; % particles size tested (m)
RhoP_test = 850:10:950; % plage de densite a tester (kg.m⁻³)

type_dict = containers.Map({'fibre','fragment','film','mousse'},0:3); % possible types
type_name = 'fibre'; % choose a type
type = type_dict(type_name); % corresponding int

nPart = 50e3; % number of particles

tf = 1e5; % simulation time (s)
dt_test = 60*60; % test time interval (s)

clear Resultats

Resultats(1:length(RhoP_test)) = struct('RhoP', 0, 'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'MeanErreur', 0);

iRes = 0;
for RhoP = RhoP_test
    iRes = iRes+1;
    [~, ~, PartPos] = varMP_model(SizePart, RhoP, type, nPart, tf, dt_test, wind, month, Lon0, Lat0, path);
    
    % 2) calculer l'erreur avec les mesures :
    % =======================================
    hcalc = histogram(PartPos, "BinEdges",  boundTest, 'Visible', 'off').Values;
    NpartModel = zeros(size(ConcentrationSample));
    j = 0;
    for i = 1:length(hcalc)
        if mod(i,2) ~= 0
            j = j+1;
            NpartModel(j) = hcalc(i);
        end
    end
    ConcentrationModel = NpartModel/dh * (nPart/L) ; % on ramène le modèle à 1
%     alpho = mean(ConcentrationSample(ConcentrationModel>0)./ConcentrationModel(ConcentrationModel>0));
%     alpha = ConcentrationSample\ConcentrationModel;
    alpha = ConcentrationModel\ConcentrationSample;
    
    
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample)./ConcentrationSample;
    meanErreur = mean(Erreur);
    
    Resultats(iRes).RhoP = RhoP;
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).Erreur = Erreur;
    Resultats(iRes).MeanErreur = meanErreur;
end

% figure(1)
% plot([Resultats.RhoP],[Resultats.Alpha])
figure(1), clf,
plot(ConcentrationSample,-Ztest_,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
xlim([0 max(ConcentrationSample)])
ylim([-L 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on
for res = Resultats
    plot(res.ConcentrationModel*res.Alpha,-Ztest_,'DisplayName', ['RhoP = ' num2str(res.RhoP)])
end
legend('Location', 'southeast')
hold off

figure(2)
hold on
plotErr = zeros(size([Resultats.RhoP]));
for i = 1:length(Resultats)
    plotErr(i) = mean(Resultats(i).Erreur);
end
plot([Resultats.RhoP], plotErr)
hold off
