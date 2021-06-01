% clear

tf = 60*60*24*5; % Simulation time
dt_test = 60*60*2;

%% Water column parameters
L = 50;
N = 50;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

%% Diffusivity
% DensiteFevrierRhoma
% KZ_day = KZ_Fev10;
% Row_day = Row_Fev10;
% z_day = z_Fev10;
% z__day = z__Fev10;
% [K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day);
% rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 
wind = 2;
date = datetime(1999,12,19);
[KZ_day,Row_day,z_day,z__day] = KsSalTemp(wind, date);
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day');
rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 

clearvars -except tf dt_test L N dz z z_ K dK rhow mean1 std1 mean2 std2

%% Particules initialisation
nPart = 10e3; % number of particles
% zPart = sparse(linspace(0, L, nPart)); % depth of particles
zPart = 20*ones(1,nPart);

sizePart1 = ones(size(zPart))*350e-6;
pd = makedist('Normal', 'mu', 350e-6, 'sigma', 50e-6);
sizePart2 = random(pd,size(zPart));

rhop = 1025;

% create list of particles
mp1 = getMPlist(nPart, sizePart1, rhop, rhow, 0);
mp2 = getMPlist(nPart, sizePart2, rhop, rhow, 0);


[zFinal1,~] = MP_simulator(mp1, zPart, K, dK, L, dz, tf, dt_test, 60*30);
[zFinal2,~] = MP_simulator(mp2, zPart, K, dK, L, dz, tf, dt_test, 60*30);

[meanConc1, stdConc1] = getMeanConc(zFinal1, length(z_), dz);
[meanConc2, stdConc2] = getMeanConc(zFinal2, length(z_), dz);

bTop = 1;
bBottom = 8;
[meanSizeT, stdSizeT,  s_, s] = getDomainSizeRep(bTop, bBottom, mp2, zFinal2);
% [meanSizeB, stdSizeB, s_, s] = getDomainSizeRep(bBottom, L, mp2, zFinal2);

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
sInit = histogram(sizePart2, 'BinEdges', s).Values;
xlabel('Particle size (µm)')
ylabel('Occurence')
title('Particle size distribution')

figure(3), clf,
hold on
plot(s_, sInit/sum(sInit), 'DisplayName', 'Initial')
plot(s_, meanSizeT/sum(meanSizeT), 'DisplayName', [num2str(bBottom) '-' num2str(bTop) ' m'] )
% plot(s_, meanSizeB/sum(meanSizeB), 'DisplayName', ['Bottom : ' num2str(bBottom) '-' num2str(L) ' m'])
hold off
legend('Location', 'best')
xlabel('Size (µm)')
xlabel('Particle size (µm)')
ylabel('Occurence')
title('Particle size distribution')

%% Kolmogorov–Smirnov test
% H0 : les répartition testées sont tirées de la même loi de distribution
% alpha = 0.05 (default)

uT = repelem(s_,round(meanSizeT));
% uB = repelem(s_,round(meanSizeB));
[h,p,ksstat] = kstest(sizePart2,'CDF',pd);
[hT,pT,ks2statT] = kstest2(sizePart2,uT);
% [hB,pB,ks2statB] = kstest2(sizePart2,uB);

[hT2,pT2,ks2statT2] = kstest(uT,'CDF',pd);

figure(4), clf,
hold on,
cdfplot(sizePart2)
cdfplot(uT)
% cdfplot(uB)
hold off,
legend('Initial',[num2str(bBottom) '-' num2str(bTop) ' m'],'Location','best')


