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
addpath('./fig-Visser/');
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day);

rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 

clearvars -except tf dt_test L N dz z z_ K dK rhow mean1 std1 mean2 std2

%% Particules initialisation
nPart = 5e3; % number of particles
zPart = linspace(0, L, nPart); % depth of particles
% sizePart = linspace(300e-6, 400e-6, nPart); % size of particle
sizePart = ones(size(zPart))*350e-6;
rhop = 1020;

% allocate memory to store particles
mp(nPart) = MP; % array of MP objects
% Fill the array
for i = 1:(nPart)
    mp(i) = MP(sizePart(i), rhop, rhow);
end, clear i,

[zFinal] = MP_simulator(mp, zPart, K, dK, L, dz, tf, dt_test, 60*30);
% zFinal = reshape(zFinal, [numel(zFinal) 1]); 

hConc = NaN(size(zFinal,1),length(z_));
for hStep = 1:size(zFinal,1)
    pp = zFinal(hStep,:);
    [histi,~] =  groupcounts(pp',z,'IncludeEmptyGroups',true);
    hConc(hStep,:) = histi'/dz*L/nPart;
end, clear hStep pp histi,
meanConc = mean(hConc, 'omitnan');
stdConc = std(hConc, 'omitnan');

plot(meanConc, -z_, meanConc+2*stdConc, -z_, '--', meanConc-2*stdConc, -z_, '--');

% [groupMP,~] =  groupcounts(zFinal,z,'IncludeEmptyGroups',true);
% conc = groupMP'/dz*L/nPart;
% plot(conc,-z_);

