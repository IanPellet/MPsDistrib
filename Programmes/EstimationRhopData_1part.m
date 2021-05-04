function [Resultats] = EstimationRhopData_1part(SizePart, type_name, RhoP_test)

%% Set Parameters

path = '../Results/EstimRho/';

wind = 1; % km/h
month = 03; 

if nargin == 0
    SizePart = false; % particles size tested (m)
    type_name = false;
    RhoP_test = 800:100:1000; % plage de densite a tester (kg.m⁻³)
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
L = 58;

[ConcentrationSample, DepthSample] = getDataNpart(type_name, SizePart, true);

dh = 0.15; % Net oppening (m)

% Boundaries in between which modeled part number have to be tested
% boundTest = zeros(1,length(DepthSample)*2);
% for i = 1:length(boundTest)
%     temp_depth = DepthSample(fix((i+1)/2));
%     if i>1 && boundTest(i-1) == temp_depth
%         boundTest(i) = temp_depth+dh;
%     else
%         boundTest(i) = temp_depth;
%     end
% end
boundTest = 0:dh:L;


Ztest_ = DepthSample + dh/2;
zplot = boundTest(1:end-1)+dh/2;

nPart = 50e3; % number of particles

tf = 1e5; % simulation time (s)
dt_test = 60*60; % test time interval (s)

clear Resultats

Resultats(1:length(RhoP_test)) = struct('RhoP', 0, 'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'MeanErreur', 0);
if SizePart
    modSize = SizePart;
else
    modSize = 1010.5e-6;
end


iRes = 0;
for RhoP = RhoP_test
    iRes = iRes+1;
    [~, ~, PartPos] = varMP_model(modSize, RhoP, TypePart, nPart, tf, dt_test, wind, month, Lon0, Lat0, path);
    
    figure(2)
    hist = histogram(PartPos,'Visible', 'off').Values;
    plot(hist,-linspace(0,L,length(hist)))
    
    % 2) calculer l'erreur avec les mesures :
    % =======================================
    hcalc = histogram(PartPos, "BinEdges",  boundTest, 'Visible', 'off').Values;
    NpartModel = hcalc(fix(DepthSample./dh+1));
%     NpartModel = zeros(size(ConcentrationSample));
%     j = 0;
%     for i = 1:length(hcalc)
%         if mod(i,2) ~= 0
%             j = j+1;
%             NpartModel(j) = hcalc(i);
%         end
%     end
    conc = hcalc/dh * (nPart/L) ;
    
    ConcentrationModel = NpartModel'/dh * (nPart/L) ; % on ramène le modèle à 1
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
    Resultats(iRes).conc = conc;
end

% typeSTR = {'fibre','fragment','film','mousse'}; % List of existing type
if type_name
    dispType = type_name;
else
    dispType = 'all';
end
ttl = ['Particle size : ' num2str(modSize*1e6) 'µm -- Particle type : ' dispType];

figure(1), clf,
plot(ConcentrationSample,-Ztest_,'pm','MarkerSize', 10, 'DisplayName', 'Sampled Data');
% xlim([0 max(ConcentrationSample)])
ylim([-L+0.75 0])
xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')
hold on
for res = Resultats
    plot(res.conc*res.Alpha,-zplot,'DisplayName', ['RhoP = ' num2str(res.RhoP)])
end
legend('Location', 'southeast')
title(ttl)
hold off

figure(2), clf,
hold on
plotErr = zeros(size([Resultats.RhoP]));
for i = 1:length(Resultats)
    plotErr(i) = mean(Resultats(i).Erreur*100);
end
plot([Resultats.RhoP], plotErr)
xlabel('Rho_p (kg.m⁻³)');
ylabel('Error (%)');
title(ttl)
hold off

end
