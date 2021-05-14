function [Resultats,ConcentrationSample,Ztest_,z_] = Estimation1Rhop_Eul(SizePart, RhoP_test, saveFig, N)
%%ESTIMATIONRHOPDATA_1PART_EUL Finds the modeled particle density closest
%%to data
% SizePart : (double) size of the modeled particles (m)
% RhoP_tes : (double array) modeled particle densities tested (kg.m⁻³)
% saveFig : (boolean) if true figures will be saved to '../Results/EstimRho/'
% N : (int) default = 500, number of meshes. High : +precise -fast / Low : -precise +fast
% 
% example : EstimationRhopData_1part_Eul(400e-6, 900:1100, true, 200)

%% Set Parameters
path = '../Results/EstimRho/'; % saved figures directory

%% Water column parameters
L = 51; % depth (m)
if nargin < 4
    N = 500; % number of meshes. High : +precise -fast / Low : -precise +fast
end
dz= L/N; % size of meshes
z = (0:dz:L)'; % meshes boundaries
z_ = z(1:end-1)+dz/2; % center of the meshes

day = '10fev'; % day corresponding to diffusive turbulence data

%% Data
% SamplingDate = datetime('3/18/2021');
% DataFile = '../Data/data_mps.txt';
% [ConcentrationSample, DepthSample] = getDataNpart(false, SizePart, true, DataFile, SamplingDate);
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

% launch simulation for each rhop
for iRes = 1:length(RhoP_test)
    RhoP = RhoP_test(iRes);    
    
    [conc, z_] = Transport_Eulerian(modSize, RhoP, N, L, day);
    
    % find the modeled concentration at depth corresponding with sample
    ConcentrationModel = zeros(size(ConcentrationSample));
    j = 0;
    for ic = 1:2:length(ibound)
        j = j+1;
        ConcentrationModel(j) = mean(conc(ibound(ic):ibound(ic+1)));
    end, clear ic j,
    
    % Compute the proportionnality coefficient minimizing the square error
    % between model and data
    alpha = ConcentrationModel\ConcentrationSample;
    
    % Compute the root mean square error between model and data
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
    
    % Fill the result structure 
    Resultats(iRes).RhoP = RhoP;
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).Erreur = Erreur;
    Resultats(iRes).rmseErreur = rmseErreur;
    Resultats(iRes).conc = conc';
end, clear iRes,

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
legend('Location', 'southeast')
title(ttl)
hold off

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
    
    fileName = ['size' num2str(modSize*1e6) '_rho' rhoInter 'part_10fev'];
    
    F = [f1, f2, f3];
    xPart = '1part_';
    N = {'profils_', 'error_', 'min_'};
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
