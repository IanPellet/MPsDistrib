clear

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

clearvars -except tf dt_test L N dz z z_ K dK rhow

%% Particules initialisation
nPart = 1e3; % number of particles
zPart = linspace(0, L, nPart); % depth of particles
sizePart = linspace(250e-6, 500e-6, nPart); % size of particle
% sizePart = ones(size(zPart))*400e-6;
rhop = 1020;

% allocate memory to store particles
mp(nPart) = MP; % array of MP objects
% Fill the array
for i = 1:(nPart)
    mp(i) = MP(sizePart(i), rhop, rhow);
end, clear i,

[zFinal] = MP_simulator(mp, zPart, K, dK, L, dz, tf, dt_test);

[groupMP,~] =  groupcounts(zFinal',z,'IncludeEmptyGroups',true);
conc = groupMP'/dz*L/nPart;
plot(conc,-z_);

