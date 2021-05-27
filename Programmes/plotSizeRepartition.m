function [resE] = plotSizeRepartition(D_test, RhoP, nDepth, saveFig)
%%PLOTLAGREUL_VARSIZE Plot size histogram in each mesh
% D_test : (double array) default = [125;225;375;750;1250;1750;2750;3750;4750]*1e-6, array of particles sizes (m)
% RhoP : (double) default = 1030, density of particles (kg.m⁻³)
% N : (int)  default = 10, number of meshes (depth of the water column = 50m)
% saveFig : (boolean) default = false, if true figures will be saved to '../Results/SizeRep/'

path = '../Results/SizeRep/';

if nargin < 4
    saveFig = false;
end
if nargin < 3
    nDepth = 10; % number of meshes
end
if nargin < 2
    RhoP = 1030; % particles density
end
if nargin < 1
    D_test = [125;225;375;750;1250;1750;2750;3750;4750]*1e-6;
end

if size(D_test,2) ~= 1
    D_test = D_test';
end

N = 200;
L = 50;
dz = L/N;

%%
dzBarPlot = L/nDepth;
zBarPlot = 0:dzBarPlot:L;
z_BarPlot = (zBarPlot(1:end-1)+zBarPlot(2:end))/2; 

%% Find corresponding depths to get from the model
ibound = NaN(size(zBarPlot)); % index boundaries of meshes to avg

for i = 1:length(ibound)
    ibound(i) = min(N,max(1,fix(zBarPlot(i)/dz)+1));
end, clear i,

%% Run simulations

resE = struct('C', [], 'z', [], 'D', 0);
for i=1:length(D_test)
    D = D_test(i);
    [resE(i).C, resE(i).z] = Transport_Eulerian(D, RhoP, N, L, '10fev');
    resE(i).D = D;
    resE(i).C = resE(i).C';
    % find the modeled concentration at depth inputed
    CModel = NaN(nDepth,1);
    j = 0;
    for ic = 1:length(ibound)-1
        j = j+1;
        CModel(j) = mean(resE(i).C(ibound(ic):ibound(ic+1)));
    end, clear ic j,
    resE(i).cBarPlot = CModel;
end, clear i,



%% Plot results
sumMPz = sum([resE.cBarPlot],2);

f1 = figure(1); clf,
barh(-z_BarPlot, [resE.cBarPlot]./sumMPz*100)
legend([num2str([resE.D]'*1e6) repmat(' µm',length(D_test),1)],'Location', 'best')
plots = get(gca, 'Children');
legend(plots)
% ylim([-5 0])
xlabel('Size repartition (%)')
ylabel('Depth (m)')
title(['Rho = ' num2str(RhoP)...
    'kg.m⁻³ -- N = ' num2str(N)] )


f2 = figure(2);clf,
hold on 
for i=1:length(D_test)
    plot([resE(i).C], -resE(i).z, 'DisplayName', ['Eul. D = ' num2str(resE(i).D*1e6) ' µm']);
end, clear i,
legend('Location', 'best')
plots = get(gca, 'Children');
legend(plots)
xlabel('Concentration (mps.m⁻³)')
ylabel('Depth (m)')
xlim([0 0.5])
ylim([-50 0])
title(['Rho = ' num2str(RhoP)...
    'kg.m⁻³ -- N = ' num2str(N)] )
hold off

%% Save figures
if saveFig
    if length(D_test) == 1
        sInter = num2str(D_test);
    else
        sInter = [num2str(min(D_test)) '-' num2str(D_test(2)-D_test(1)) '-' num2str(max(D_test))];
    end
    
    fileName = ['size' sInter '_rho' num2str(RhoP)];
    
    F = [f1, f2];
    nFiles = {'rep_', 'profiles_'};
    for i=1:length(F)
        f = F(i);
        n = nFiles{i};
        NamePNG = [path 'png/' n fileName '.png'];
        NameFIG = [path 'fig/' n fileName '.fig'];
        exportgraphics(f,NamePNG);
        savefig(f,NameFIG);
    end, clear i,
end

end