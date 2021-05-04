function [Resultats] = EstimationRhopData_2part(SizePart, type_name, RhoP_test)

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
L = 55.1;

[ConcentrationSample, DepthSample] = getDataNpart(TypePart, SizePart, true);

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


if SizePart
    modSize = SizePart;
else
    modSize = 1010.5e-6;
end


%% Run model
clear Resultats

Resultats(1:length(RhoP_test)) = struct('Rho1', 0, 'Rho2', 0, 'PartPos',...
    [] ,'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'MeanErreur', 0);


combRho = nchoosek(RhoP_test,2);

for iRes=1:size(combRho,1)
    Rho1 = RhoP_test(iRes,1);
    Rho2 = RhoP_test(iRes,2);
    
    [~, ~, PartPos1] = varMP_model(modSize, Rho1, TypePart, nPart, tf, dt_test, wind, month, Lon0, Lat0, path);
    [~, ~, PartPos2] = varMP_model(modSize, Rho2, TypePart, nPart, tf, dt_test, wind, month, Lon0, Lat0, path);
    
    Resultats(iRes).Rho1 = Rho1;
    Resultats(iRes).Rho2 = Rho2;
    Resultats(iRes).PartPos = [PartPos1 PartPos2];
end, clear iRes,

%% Model output treatment
for i=1:length(Resultats)
    hcalc = histogram(Resultats(i).PartPos, "BinEdges",  boundTest, 'Visible', 'off').Values;
    NpartModel = hcalc(fix(DepthSample./dh+1));
    conc = hcalc/dh * (nPart/L) ;
    ConcentrationModel = NpartModel'/dh * (nPart/L) ; % on ramène le modèle à 1

    alpha = ConcentrationModel\ConcentrationSample;
   
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample)./ConcentrationSample;
    meanErreur = mean(Erreur);
    
    Resultats(i).ConcentrationModel = ConcentrationModel;
    Resultats(i).Alpha = alpha;
    Resultats(i).Erreur = Erreur;
    Resultats(i).MeanErreur = meanErreur;
    Resultats(i).conc = conc;
end, clear i,


%% Display results
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
    plot(res.conc*res.Alpha,-zplot,'DisplayName', ['Rho1 = ' num2str(res.Rho1) ', Rho2 = ' num2str(res.Rho2)])
end
legend('Location', 'southeast')
title(ttl)
hold off

% figure(2), clf,
% hold on
% plotErr = zeros(size([Resultats.RhoP]));
% for i = 1:length(Resultats)
%     plotErr(i) = mean(Resultats(i).Erreur*100);
% end
% plot([Resultats.RhoP], plotErr)
% xlabel('Rho_p (kg.m⁻³)');
% ylabel('Error (%)');
% title(ttl)
% hold off

end
