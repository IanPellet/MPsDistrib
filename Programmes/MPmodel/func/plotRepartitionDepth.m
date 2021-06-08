IDref = '2021-6-4-14-47-29_686516';
ID = '2021-6-7-15-7-28_342636';

load(['/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/' IDref '.mat'], 'meanConc')
meanConcRef = cell2mat(meanConc);
load(['/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/' ID '.mat'])
meanConc = cell2mat(meanConc);


dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  
nProfile = length(meanConc);

[~,Row_day,~,z__day] = KsSalTemp(wind, date);
rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 
mp = getMPlist(nPart, sizeP, rhop, rhow, 0);

f1 = figure(1); clf, hold on,
plot(meanConcRef(end,:), -z_, '--', 'DisplayName', 'Constant Size Repartition')
plot(meanConc(end,:), -z_, 'DisplayName', 'Uniform Size Repartition')
hold off
legend('Location', 'best')
xlabel('Concentration (mps.m⁻³)')
ylabel('Depth (m)')

% deltaConc = abs(meanConc-meanConcRef)./meanConcRef;
deltaConc = (meanConc-meanConcRef);
deltaConc(isnan(deltaConc)) = 0;
% Plot deltaConc after t = 100h
tStab = 100*60*60;
iStab = tStab/dtAvgC;

figure(2), clf,
plot(deltaConc(iStab:end,:)*100, -z_)
xlabel('deltaConcentration (mps.m⁻³)')
ylabel('Depth (m)')
title('Average concentration differance between simulation with constant and variable particle size',...
    ['Computed each ' num2str(dtAvgC/60) ' min since t_s_t_a_b = ' num2str(tStab/60/60) ' h'])

meanDeltaConc = mean(deltaConc(iStab:end,:));
stdDeltaConc = std(deltaConc(iStab:end,:),0);
f3 = figure(3); clf,
hold on
plot(meanDeltaConc*100, -z_)
plot((meanDeltaConc+2*stdDeltaConc)*100, -z_, '--')
plot((meanDeltaConc-2*stdDeltaConc)*100, -z_, '--')
xlabel('deltaConcentration (mps.m⁻³)')
ylabel('Depth (m)')
title('Mean concentration profile difference +/- 2 std', ['from t_s_t_a_b = '...
    num2str(tStab/60/60) ' h to t_f_i_n_a_l = ' num2str(tf/60/60) ' h'])


%% Compute time average of size repartition
load(['/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/Data/' ID '-zHist864000.mat'])


% bTopT = 0;
% bBottomT = 10;
% bTopB = 60;
% bBottomB = L;
% [meanSizeT, stdSizeT, meanNT, ~] = getDomainSizeCDF(bTopT, bBottomT, mp, zHistory(end-60*60/dt:end));
% [meanSizeB, stdSizeB, meanNB, s] = getDomainSizeCDF(bTopB, bBottomB, mp, zHistory(end-60*60/dt:end));
% [meanSizeM, stdSizeM, meanNM, ~] = getDomainSizeCDF(bBottomT, bTopB, mp, zHistory(end-60*60/dt:end));

bound = 0:5:L;
bound_ = (bound(1:end-1) + bound(2:end))/2;
[meanSize1, stdSize1, meanN1, s] = getDomainSizeCDF(bound(1), bound(2), mp, zHistory(end-60*60/dt:end));
meanSize = NaN(length(bound)-1,length(meanSize1));
stdSize = NaN(length(bound)-1,length(stdSize1));
meanN = NaN(length(bound)-1,1);
meanSize(1,:) = meanSize1;
stdSize(1,:) = stdSize1;
meanN(1) = meanN1;
for i = 2:length(bound)-1
    [meanSize(i,:), stdSize(i,:), meanN(i), ~] = getDomainSizeCDF(bound(i), bound(i+1), mp, zHistory(end-60*60/dt:end));
end, clear i,

% clear zHistory mp,


[iecdf,x] = ecdf(sizeP);
[C,ia,~] = unique(x);
sizeCDFinit =interp1(C,iecdf(ia),s);

clear x iecdf C ia 

f4 = figure(4); clf,
hold on
plot(s,sizeCDFinit, 'DisplayName', 'Initial size CDF')
for i = 1:length(meanN)
    plot(s,meanSize(i,:), 'DisplayName', ['Size CDF : ' num2str(bound(i)) '-' num2str(bound(i+1)) ' m'])  
end, clear i,
hold off
legend('Location', 'best')
xlabel('Size (m)')
ylabel('Cumulative probability')
lines = get(gca, 'Children');
set(lines, {'Color'}, nToColorMap(length(lines)))

%% Kolmogorov–Smirnov test
rej = NaN(length(bound_),1);
D = NaN(length(bound_),1);
for i = 1:length(bound_)
    [rej(i), D(i)] = MP_kstest(sizeCDFinit, nPart, meanSize(i,:), meanN(i), 0.05);
end, clear i,

%% Save fig
figName = [path ID '-VS-' IDref '-C.fig'];
savefig(f1, figName);
figName = [path ID '-VS-' IDref '-deltaC.fig'];
savefig(f3, figName);
figName = [path ID '-sizeRep.fig'];
savefig(f4, figName);
        

