function [resE] = plotLagrEul_varSize(D_test, RhoP, N)
%%PLOTLAGREUL_VARSIZE Plot size histogram in each mesh
% D_test : array of particles sizes (m)
% RhoP : density of particles
% N : number of meshes (depth of the water column = 50m)
    
if nargin == 0
    %% Choose the MPs sizes : D_test

    % DataFile = '../Data/data_mps.txt';
    % dataTable = load_MPs_data(DataFile);
    % D_test = unique(dataTable(:,'med_size').Variables)*1e-6; % sorted list of depth

    D_test = [125
    225
    375
    750
    1250
    1750
    2750
    3750
    4750
    ]*1e-6;

%     D_test = linspace(100e-6,5e-3,10)';

    % D_test = 350e-6;

    %% Parameters 
    RhoP = 1010.05; % particles density
    N = 50; % number of meshes (water column depth = 50)
end

if size(D_test,2) ~= 1
    D_test = D_test';
end

%% Run simulations
resE = struct('C', [], 'z', [], 'D', 0);
for i=1:length(D_test)
    D = D_test(i);
    [resE(i).C, resE(i).z] = Transport_Eulerian(D, RhoP, N);
    resE(i).D = D;
    resE(i).C = resE(i).C';
end, clear i,

%% Plot results
sumMPz = sum([resE.C],2);

figure(1); clf,
barh(-resE(1).z, [resE.C]./sumMPz*100)
legend([num2str(D_test*1e6) repmat(' µm',length(D_test),1)],'Location', 'best')
% ylim([-5 0])
xlabel('Size repartition (%)')
ylabel('Depth (m)')


figure(2);clf,
hold on 
for i=1:length(D_test)
%     alpha = (sum(resL(i).C)*dz_L)/(sum(resE(i).C)*dz_E);
    plot([resE(i).C], -resE(i).z, 'DisplayName', ['Eul. D = ' num2str(resE(i).D*1e6) ' µm']);
%     plot([resL(i).C], -resL(i).z, 'DisplayName', ['Lag. D = ' num2str(resL(i).D*1e6) ' µm']);
%     plot([resE(i).C]*alpha, -resE(i).z, 'DisplayName', ['Eulerian model']);
%     plot([resL(i).C], -resL(i).z, 'DisplayName', ['Lagrangian model']);
end, clear i,
legend('Location', 'best')
xlabel('Concentration (mps.m⁻³)')
ylabel('Depth (m)')
xlim([0 0.5])
ylim([-50 0])
% title(['D = ' num2str(D_test(1)*1e6) 'µm -- rho = ' num2str(rhoP)...
%     'kg.m⁻³ -- nPart = ' num2str(nPart) ' -- N = ' num2str(N)] )
title(['Rho = ' num2str(RhoP)...
    'kg.m⁻³ -- N = ' num2str(N)] )
hold off


end