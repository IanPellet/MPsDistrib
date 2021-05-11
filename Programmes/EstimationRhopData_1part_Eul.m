function [Resultats,ConcentrationSample,Ztest_,z_] = EstimationRhopData_1part_Eul(SizePart, RhoP_test, saveFig)
%%ESTIMATIONRHOPDATA_1PART_EUL Finds the modeled particle density closest
%%to data
% SizePart : (double) size of the modeled particles (m)
% RhoP_tes : (double array) modeled particle densities tested (kg.m⁻³)
% saveFig : (boolean) if true figures will be saved to '../Results/EstimRho/'


%% Set Parameters
path = '../Results/EstimRho/'; % saved figures directory

% Wind parameters
wind = 50; % km/h
month = 03; 
day = '10fev';


%% Data loading
L = 55;
N = L;
dz= L/N;
z = 0:dz:L;
z_ = z(1:end-1)+dz/2;

% SamplingDate = datetime('3/18/2021');
% DataFile = '../Data/data_mps.txt';
% [ConcentrationSample, DepthSample] = getDataNpart(false, SizePart, true, DataFile, SamplingDate);
CMes=[0.62 0.34 0.06 0.02 0];
ZMes=[1 10 15 40 L];
ConcentrationSample = interp1(ZMes,CMes,z_,'pchip')';
DepthSample = z_';

dh = 0.71; % Net oppening (m)

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

Ztest_ = DepthSample - dh/2;

clear Resultats

Resultats(1:length(RhoP_test)) = struct('RhoP', 0, 'ConcentrationModel', zeros(size(ConcentrationSample)),...
    'Alpha', 0, 'Erreur', zeros(size(ConcentrationSample)), 'rmseErreur', 0);

if SizePart
    modSize = SizePart;
else
    modSize = 350e-6;
end

for iRes = 1:length(RhoP_test)
    RhoP = RhoP_test(iRes);    
    [conc, z_] = Transport_Eulerian(modSize, RhoP, N, day);
    
    ConcentrationModel = conc';

    alpha = ConcentrationModel\ConcentrationSample;
    
    
    Erreur = abs(ConcentrationModel.*alpha - ConcentrationSample);
    rmseErreur = sqrt(mean(Erreur.^2,'omitnan'));
    
    Resultats(iRes).RhoP = RhoP;
    Resultats(iRes).ConcentrationModel = ConcentrationModel;
    Resultats(iRes).Alpha = alpha;
    Resultats(iRes).Erreur = Erreur;
    Resultats(iRes).rmseErreur = rmseErreur;
    Resultats(iRes).conc = conc;
end, clear iRes,


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
        NameFIG = [path 'fig/' xPart n fileName '.fig'];
        exportgraphics(f,NamePNG);
        savefig(f,NameFIG);
    end, clear i,
end

end
