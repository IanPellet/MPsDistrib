chemin = '../Results/frag/';
runID = strrep(join(split(num2str(clock)),'-'),'.','_'); % unique ID to name files
runID = runID{:};

pc = 0.490; % plane of weakness condition
D = log(8*pc)/log(8);
% r = 1e-4:1e-6:5e-3;
% r = 200e-6:10e-6:10e-2;
r = logspace(-4,-1,1e2);
r_ = (r(1:end-1) + r(2:end))/2;
Nn = r.^(-D);
nPartTarg = 50e3;
Nn = Nn./sum(Nn)*nPartTarg;
Nn_ = (Nn(1:end-1) + Nn(2:end))/2;

sizeP = repelem(r,round(Nn));

sN = histogram(sizeP, 'BinEdges', r).Values;
hold on
plot(r,Nn)
hold off

loglog(r_, sN./diff(r),'+')
hold on
loglog(r_, Nn_./diff(r),'+')
hold off

nPart = length(sizeP);

%% Time parameters
dt_test = 60*60*2;
tf = 60*60*24*3;

%% Water column parameters
%% Load hydrodinamic model data
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro, 'H0')

%% Water column parameters
% Find depth of the column
Station = 'RN2';
% Load file with indexes corresponding to each stations
stationFile = '../Data/stationIJ_CEREGE.mat';
load(stationFile,'stationIJ');
% Get the indexes of the right station
I0 = stationIJ{stationIJ{:,'station'} == Station,'I0'};
J0 = stationIJ{stationIJ{:,'station'} == Station,'J0'};
L = H0(I0,J0); % depth
clear H0  stationIJ,

N = 50;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh 

zInit = random('Uniform', 0,L,1,nPart);
frag = 0;

%% Diffusivity
date = datetime(2020,03,18);
wind = 30;
rhop = 1025;
[KZ_day,Row_day,z_day,z__day] = KsSalTemp(wind, date);
[K,dK] = Diffusivity(z,z_,dz,0.8,0,KZ_day,z_day');
rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 
clear KZ_day Row_day z_day z__day,


%% Create particles
disp('Create Particles')
mp = getMPlist(nPart, sizeP, rhop, rhow, frag);


%% Simulation
dtAvgC = 30*60;

disp('Run Simulator')
[zFinal, dt] = MP_simulator(mp, zInit, K, dK, L, dz, tf, dt_test, dtAvgC);

disp('Save parameters')
save([chemin runID '_param.mat'],...
    'runID', 'dt_test', 'date', 'tf', 'L', 'N', 'dz', 'nPart',...
    'sizeP', 'zInit', 'frag', 'wind', 'rhop', 'chemin', 'r',...
    'K', 'dK', 'dtAvgC', 'rhow', 'mp', 'pc');

clearvars -except zFinal mp dt L runID chemin,

disp('Save results')
save([chemin runID '_res.mat'], 'dt', 'zFinal');


%% Compute time average of size repartition (CDF)
bound = 0:5:L;
% bound_ = (bound(1:end-1) + bound(2:end))/2;
clear L dt,
% [meanSize1, stdSize1, meanN1, s] = getDomainSizeCDF(bound(1), bound(2), mp, zFinal);
% 
% parameters and 
% meanSizeCDF = NaN(length(bound)-1,length(meanSize1));
% stdSizeCDF = NaN(length(bound)-1,length(stdSize1));
% meanN = NaN(length(bound)-1,1);
% meanSizeCDF(1,:) = meanSize1;
% stdSizeCDF(1,:) = stdSize1;
% meanN(1) = meanN1;
% for i = 2:length(bound)-1
%     [meanSizeCDF(i,:), stdSizeCDF(i,:), meanN(i), ~] = getDomainSizeCDF(bound(i), bound(i+1), mp, zFinal);
% end, clear i,
% meanSizeCDF = meanSizeCDF(sum(meanSizeCDF,2, 'omitnan')~=0,:);
% % clear zHistory mp,
% 
% 
% [iecdf,x] = ecdf(sizeP);
% [C,ia,~] = unique(x);
% sizeCDFinit =interp1(C,iecdf(ia),s);
% 
% clear x iecdf C ia 
% 
% f4 = figure(4); clf,
% hold on
% plot(s,sizeCDFinit, 'DisplayName', 'Initial size CDF')
% for i = 1:size(meanSizeCDF,1)
%     plot(s,meanSizeCDF(i,:), 'DisplayName', ['Size CDF : ' num2str(bound(i)) '-' num2str(bound(i+1)) ' m'])  
% end, clear i,
% hold off
% legend('Location', 'best')
% xlabel('Size (m)')
% ylabel('Cumulative probability')
% lines = get(gca, 'Children');
% set(lines, {'Color'}, nToColorMap(length(lines)))

%% Size rep // power law
topBound = 0; bottomBound = 5;
% Get indexes of particles in the domain
zFmat = cell2mat(zFinal);
% clear zFinal,
iDomain = zFmat >= topBound & zFmat <= bottomBound;
nPartDom = sum(iDomain,2); % number of particles in the domain at each time step
avgNDom = mean(nPartDom); % time avg of number of part in the domain

% sizeCat = logspace(r(1),r(end),50);
sizeCat = logspace(-4,-1,50);

nCat = length(sizeCat);
histi = NaN(size(zFmat,1),nCat-1);
for iz=1:size(zFmat,1)

%     ppSort = [mp(iDomain(iz,:)).size_];
% 
%     topBound = sizeCat(2);
%     j = 1;
%     count = zeros(nCat,1);
% 
%     for i=1:length(ppSort)
%         if ppSort(i) > topBound
%             
%             topBound = sizeCat(j+2);
%             j = j+1;
%             
%         end
%         count(j) = count(j)+1;
%     end, clear i,
% 
%     histi(iz,:) = count(1:end-1);
    histi(iz,:) = histogram([mp(iDomain(iz,:)).size_], 'BinEdges', sizeCat).Values;
end, clear iz,

% ppSort = sizeP;
% 
% topBound = sizeCat(2);
% j = 1;
% count = zeros(nCat,1);
% 
% for i=1:length(ppSort)
%     if ppSort(i) > topBound
% 
%         topBound = sizeCat(j+2);
%         j = j+1;
% 
%     end
%     count(j) = count(j)+1;
% end, clear i,
% 
% histINIT = count(1:end-1);
histINIT = histogram([mp.size_], 'BinEdges', sizeCat).Values;
clear sizeP,

meanHist = mean(histi);
stdHist = std(histi);

sampDepth = ['_' num2str(topBound) '-' num2str(bottomBound) 'm'];

disp('Save parameters and results')
save([chemin runID '-sizeRep' sampDepth '.mat'],...
    'sizeCat', 'histi', 'histINIT', 'meanHist', 'stdHist', 'bound');

sizeCat_ = (sizeCat(1:end-1) + sizeCat(2:end))/2;

delta = abs(meanHist-histINIT);
diverg = delta > 2*stdHist;
pts = meanHist./diff(sizeCat);

f1 = figure(1);
clf,
loglog(sizeCat_*1e3,histINIT./diff(sizeCat), 'DisplayName', 'Initial size repartition')
hold on
sRep = loglog(sizeCat_*1e3,meanHist./diff(sizeCat), '+', 'DisplayName', ['Final size repartition ('...
    num2str(topBound) '-' num2str(bottomBound) 'm)']);
loglog(sizeCat_*1e3,(meanHist + 2*stdHist)./diff(sizeCat), '--', 'Color', sRep.Color,...
    'DisplayName', 'Final size repartition + 2std')
loglog(sizeCat_*1e3,(meanHist - 2*stdHist)./diff(sizeCat), '--', 'Color', sRep.Color,...
    'DisplayName', 'Final size repartition - 2std')
loglog(sizeCat_(diverg)*1e3,pts(diverg), '*', 'DisplayName', 'Divergent points')
xline(5, '--', 'DisplayName', '5 mm')
xline(1, '--', 'DisplayName', '1 mm')
hold off
legend('Location', 'best')
xlabel('Size (mm)')
ylabel('Normalized abundance of MPs')
xlim([0.2 100])

n_rInit = -diff(histINIT)./diff(sizeCat_);
n_r =  -diff(meanHist)./diff(sizeCat_);
f2 = figure(2);
clf,
loglog(sizeCat__,n_rInit)


savefig(f1, [chemin runID '-sizeRep' sampDepth '.fig'])

