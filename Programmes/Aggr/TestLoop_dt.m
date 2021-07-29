dtTheo_test = logspace(-6,0,50);

dt_test = 1e-7;
date = datetime(2020,03,18);


%% Water column parameters
% Find depth of the column
load('../Data/2020waterCol_RN2.mat', 'H0')
L = H0; % depth
clear H0,

N = 20;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

nMP = 1e3;
sizeMP = ones(nMP,1)*450e-6;
zInitMP = linspace(29.9,30.1,nMP);
rhoMP = 1000;

nAgg = 1;
sizeAgg = ones(nAgg,1)*1e-3;
zInitAgg = 30;
rhoAgg = 1040;

wind = 50;

[KZ_day,Row_day,z_day,z__day] = KsSalTemp2020(wind, date);
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day');
K = 0.01*ones(size(K));
dK = zeros(size(dK));
rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 
clear KZ_day Row_day z_day z__day,

%% create list of MP
% allocate memory to store particles
mp_list(nMP) = MP; % array of MP objects
% Fill the array
for i = 1:(nMP)
    mp_list(i) = MP(sizeMP(i), rhoMP, rhow, i);
end, clear i,

%% create list of Agrr
% allocate memory to store particles
agg_list(nAgg) = Aggr; % array of Aggr objects
% Fill the array
for i = 1:(nAgg)
        agg_list(i) = Aggr(sizeAgg(i), rhoAgg, rhow, i, 0);
end, clear i,   

path = '../Results/Aggr/';



%% RUN TEST

savedCollisions = cell(length(dtTheo_test),2);

for idt = 1:length(dtTheo_test)
    disp([num2str(idt) '/' num2str(length(dtTheo_test))])
    dtTheoB = dtTheo_test(idt);
    
    [collisionB, dtB] = Part_SimulatorBoom(mp_list, agg_list, zInitMP, zInitAgg, K, dK, L, dz, dtTheoB, dt_test, dtTheoB);
    
    if dtTheoB ~= dtB
        disp(['Error test ' num2str(idt)])
    end
    
    colTestB = cell2mat(collisionB');   
    DataColB = colTestB(:,1);
    clear colTestB collisionB,
    
    collisionA = Collision_Simulator(mp_list, agg_list, zInitMP, zInitAgg, K, dK, L, dz, dtB, dt_test, dtB);
    colTestA = cell2mat(collisionA');   
    DataColA = colTestA(:,1);
    clear colTestA collisionA,
    
    savedCollisions(idt,1) = {DataColA};
    savedCollisions(idt,2) = {DataColB};
    clear DataColB DataColA,
end, clear idt,

save([path 'TestLoop_dt' num2str(length(dtTheo_test)) '_collisions_K' num2str(mean(K)) '2.mat'])