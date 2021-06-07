IDref = '2021-6-4-14-47-29_686516';
ID = '2021-6-2-20-47-30_369244';

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
plot(meanConc(end,:), -z_, 'DisplayName', 'Normal Size Repartition')
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

bTopT = 0;
bBottomT = 5;
bTopB = 10;
bBottomB = 20;
[meanSizeT, stdSizeT, meanNT, ~] = getDomainSizeCDF(bTopT, bBottomT, mp, zHistory(end-60*60/dt:end));
[meanSizeB, stdSizeB, meanNB, s] = getDomainSizeCDF(bTopB, bBottomB, mp, zHistory(end-60*60/dt:end));
[meanSizeM, stdSizeM, meanNM, ~] = getDomainSizeCDF(bBottomT, bTopB, mp, zHistory(end-60*60/dt:end));

% clear zHistory mp,


[iecdf,x] = ecdf(sizeP);
[C,ia,~] = unique(x);
sizeCDFinit =interp1(C,iecdf(ia),s);

clear x iecdf C ia 

f4 = figure(4); clf,
hold on
plot(s,sizeCDFinit, 'DisplayName', 'Initial size CDF')
plot(s,meanSizeT, 'DisplayName', ['Size CDF : ' num2str(bTopT) '-' num2str(bBottomT) ' m'])
plot(s,meanSizeM, 'DisplayName', ['Size CDF : ' num2str(bBottomT) '-' num2str(bTopB) ' m'])
plot(s,meanSizeB, 'DisplayName', ['Size CDF : ' num2str(bTopB) '-' num2str(bBottomB) ' m'])
hold off
legend('Location', 'best')
xlabel('Size (m)')
ylabel('Cumulative probability')


%% Kolmogorov–Smirnov test
[Trej, TD] = MP_kstest(sizeCDFinit, nPart, meanSizeT, meanNT, 0.05);
[Brej, BD] = MP_kstest(sizeCDFinit, nPart, meanSizeB, meanNB, 0.05);
[Mrej, MD] = MP_kstest(sizeCDFinit, nPart, meanSizeM, meanNM, 0.05);

%% Save fig
figName = [path ID '-VS-' IDref '-C.fig'];
savefig(f1, figName);
figName = [path ID '-VS-' IDref '-deltaC.fig'];
savefig(f3, figName);
figName = [path ID '-sizeRep.fig'];
savefig(f4, figName);
        

