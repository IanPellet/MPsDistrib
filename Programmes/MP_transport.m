% clear

tf = 60*60*24*2; % Simulation time
dt_test = 60*60*2;

%% Water column parameters
L = 50;
N = 50;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

%% Diffusivity
DensiteFevrierRhoma
KZ_day = KZ_Fev10;
Row_day = Row_Fev10;
z_day = z_Fev10;
z__day = z__Fev10;
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day);

rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 

clearvars -except tf dt_test L N dz z z_ K dK rhow mean1 std1 mean2 std2

%% Particules initialisation
nPart = 50e3; % number of particles
zPart = sparse(linspace(0, L, nPart)); % depth of particles

sizePart1 = ones(size(zPart))*350e-6;
pd = makedist('Normal', 'mu', 350e-6, 'sigma', 50e-6);
sizePart2 = random(pd,size(zPart));

rhop = 1020;

% create list of particles
mp1 = getMPlist(nPart, sizePart1, rhop, rhow, 0);
mp2 = getMPlist(nPart, sizePart2, rhop, rhow, 0);


[zFinal1,~] = MP_simulator(mp1, zPart, K, dK, L, dz, tf, dt_test, 60*30);
[zFinal2,~] = MP_simulator(mp2, zPart, K, dK, L, dz, tf, dt_test, 60*30);

[meanConc1, stdConc1] = getMeanConc(zFinal1, z, z_, dz, L);
[meanConc2, stdConc2] = getMeanConc(zFinal2, z, z_, dz, L);


% Plot concentration profiles
figure(1), clf, hold on,

p1 = plot(meanConc1, -z_, 'b');
plot(meanConc1+2*stdConc1, -z_, '--b', meanConc1-2*stdConc1, -z_, '--b');
p2 = plot(meanConc2, -z_, 'r');
plot(meanConc2+2*stdConc2, -z_, '--r', meanConc2-2*stdConc2, -z_, '--r');

legend([p1, p2],{'cst size','normal distrib of sizes'}, 'Location', 'best');

xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')

hold off

% Plot total size repartition
figure(2), clf,
histogram(sizePart2*1e6);
xlabel('Particle size (µm)')
ylabel('Occurence')
title('Particle size distribution')


% function [mp_list] = getMPlist(nPart, sizePart, rhop, rhow, rFrag)
%     % allocate memory to store particles
%     mp_list(nPart) = MP; % array of MP objects
%     % Fill the array
%     for i = 1:(nPart)
%         mp_list(i) = MP(sizePart(i), rhop, rhow, rFrag);
%     end, clear i,
% end
% 
% function [meanConc, stdConc] = getMeanConc(zPart, z, z_, dz, L)
%     hConc = NaN(size(zPart,1),length(z_));
%     for hStep = 1:size(zPart,1)
%         pp = zPart{hStep};
%         [histi,~] =  groupcounts(pp',z,'IncludeEmptyGroups',true);
%         hConc(hStep,:) = histi'/dz*L/numel(pp);
%     end, clear hStep pp histi,
%     meanConc = mean(hConc, 'omitnan');
%     stdConc = std(hConc, 'omitnan');
% end

