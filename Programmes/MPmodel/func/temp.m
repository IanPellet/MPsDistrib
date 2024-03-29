IDs = {'2021-6-8-14-45-56_784857'
'2021-6-8-14-53-41_974487'
'2021-6-8-15-6-21_833235'
'2021-6-8-15-14-26_722505'
'2021-6-8-15-42-18_158014'
'2021-6-8-15-51-37_84898'
'2021-6-8-16-10-52_993057'
'2021-6-8-16-18-23_554119'
% '2021-6-9-8-23-43_8441'
% '2021-6-9-8-32-45_173065'
};
% refs = {'2021-6-8-15-6-21_833235'
% '2021-6-8-15-14-26_722505'};

for iID = 1:length(IDs)
    clearvars -except iID IDs refs,
    ID = IDs{iID};
%     ID = '2021-6-9-8-23-43_8441';
%     IDref = refs{[mod(iID,2)==0]+1};
%     load(['/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/' IDref '.mat'], 'meanConc', 'stdConc')
%     meanConcRef = cell2mat(meanConc);
%     stdConcRef = cell2mat(stdConc);
%     load(['/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/' ID '.mat'])
%     meanConc = cell2mat(meanConc);
%     stdConc = cell2mat(stdConc);
% 
%     dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
%     z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  
%     nProfile = length(meanConc);
% 
%     [~,Row_day,~,z__day] = KsSalTemp(wind, date);
%     rhow = interp1(-z__day,Row_day,z,'pchip'); % density of sea water 
%     mp = getMPlist(nPart, sizeP, rhop, rhow, 0);
% 
%     %% Stab
%    f1 = figure(1); clf, hold on,
% 
% pMCref = plot(meanConcRef(end,:), -z_, 'DisplayName', 'Reference');
% pMC = plot(meanConc(end,:), -z_, 'DisplayName', 'Tested Size Repartition');
% 
% plot(meanConcRef(end,:)+2*stdConcRef(end,:), -z_, '--', 'Color', pMCref.Color, 'DisplayName', 'Reference + 2std');
% plot(meanConcRef(end,:)-2*stdConcRef(end,:), -z_, '--', 'Color', pMCref.Color, 'DisplayName', 'Reference - 2std');
% 
% plot(meanConc(end,:)+2*stdConc(end,:), -z_, '--', 'Color', pMC.Color, 'DisplayName', 'Tested Size Repartition + 2std');
% plot(meanConc(end,:)-2*stdConc(end,:), -z_, '--', 'Color', pMC.Color, 'DisplayName', 'Tested Size Repartition - 2std');
% 
% hold off
% legend('Location', 'best')
% xlabel('Concentration (mps.m⁻³)')
% ylabel('Depth (m)')
% title("Last average concentration profile", ['Average on the last ' num2str(dtAvgC/30) ' min']) 
% 
% 
%     figName = [path ID '-VS-' IDref '-C'];
%     savefig(f1, [figName '.fig']);
%     exportgraphics(f1, [figName '.png']);

    
    load(['/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/' ID '.mat'])
   
    dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
    z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  
%     nProfile = length(meanConc);
    
%     % number of indexes corresponding to testBetween in meanConc
%     diStab = testStab/dtAvgC; 
% 
%     % Times at which we compute a value for StabC
%     tStabC = (tC(1:end-diStab)+tC(1+diStab:end))/2;
% 
% StabC = NaN(1,length(meanConc)-diStab);
% stdStab = NaN(1,length(meanConc)-diStab);
% for idC = 1:length(meanConc)-diStab
%     % Compute root mean square error
% %     StabC(idC) = sqrt(mean((meanConc{idC}-meanConc{idC+diStab}-stdConc{idC+diStab}).^2));
% %     StabC(idC) = mean(abs(meanConc{idC}-meanConc{idC+diStab}));
%     a = abs(meanConc{idC}-meanConc{idC+diStab})-stdConc{idC+diStab};
%     a(a<0)=0;
%     StabC(idC) = mean(a);
%     stdStab(idC) = std(a);
% end, clear idC,
% 
% stdConcMat = cell2mat(stdConc);
% RMSstd = sqrt(mean(stdConcMat.^2,2));
% 
% f3 = figure(3); clf,
% hold on
% plot(tStabC/60/60,StabC, 'DisplayName',...
%     ['RMS(C(t) - C(t+' num2str(testStab/60/60) 'h)+/-2std)'])
% % plot(tStabC/60/60,StabC + stdStab, '--',...
% %     'DisplayName', ['RMS(std_C(t+' num2str(testStab/60/60) 'h))'])
% hold off
% legend('Location', 'best')
% xlabel('simulation time (h)')
% ylabel('Difference between the profiles (mps.m⁻³)')
% title('Evolution of the stability of the concentration profile over time',...
%     ['Tested between profiles separated by ' num2str(testStab/60/60) ' h of simulation'])
% 
% 
% 
% idC = length(meanConc)-diStab;
% 
% % figure(2), clf, 
% % hold on
% % plot(meanConc{idC}, -z_)
% % plot(meanConc{idC+diStab}, -z_)
% 
% figure(2); clf, hold on,
% 
% pMCref = plot(meanConc{idC+diStab}, -z_, 'DisplayName', 'Constant Size Repartition');
% pMC = plot(meanConc{idC}, -z_, 'DisplayName', 'Uniform Size Repartition');
% 
% plot(meanConc{idC+diStab}+2*stdConc{idC+diStab}, -z_, '--', 'Color', pMCref.Color, 'DisplayName', 'Constant Size Repartition + 2std');
% plot(meanConc{idC+diStab}-2*stdConc{idC+diStab}, -z_, '--', 'Color', pMCref.Color, 'DisplayName', 'Constant Size Repartition - 2std');
% 
% plot(meanConc{idC}+2*stdConc{idC}, -z_, '--', 'Color', pMC.Color, 'DisplayName', 'Uniform Size Repartition + 2std');
% plot(meanConc{idC}-2*stdConc{idC}, -z_, '--', 'Color', pMC.Color, 'DisplayName', 'Uniform Size Repartition - 2std');
% 
% hold off
% legend('Location', 'best')
% xlabel('Concentration (mps.m⁻³)')
% ylabel('Depth (m)')
% title("Last average concentration profile", ['Average on the last ' num2str(dtAvgC/30) ' min']) 
% 
% 
% %     StabC = NaN(1,length(meanConc)-diStab);
% %     for idC = 1:length(meanConc)-diStab
% %         % Compute root mean square error
% %         StabC(idC) = sqrt(mean((meanConc{idC}-meanConc{idC+diStab}).^2));
% %     %     StabC(idC) = mean(abs(meanConc{idC}-meanConc{idC+diStab}));
% %     end, clear idC,
% 
%     
%     [StabC,tStabC] = getStability(meanConc, testStab, dtAvgC, tC);
%     stdConcMat = cell2mat(stdConc);
%     RMSstd = sqrt(mean(stdConcMat.^2,2));
% 
%     f1 = figure(1); clf,
%     hold on
%     plot(tStabC/60/60,StabC, 'DisplayName',...
%         ['RMS(C(t) - C(t+' num2str(testStab/60/60) 'h))'])
%     plot(tStabC/60/60,RMSstd(testStab/dtAvgC+1:end), '--',...
%         'DisplayName', ['RMS(std_C(t+' num2str(testStab/60/60) 'h))'])
%     hold off
%     legend('Location', 'best')
%     xlabel('simulation time (h)')
%     ylabel('Root Mean Square Difference (mps.m⁻³)')
%     title('Evolution of the stability of the concentration profile over time',...
%         ['Tested between profiles separated by ' num2str(testStab/60/60) ' h of simulation'])
    

concMat = cell2mat(meanConc);
meanConTab = array2table(concMat, 'VariableNames', split(num2str(z_)),...
    'RowNames', split(num2str(tC)));
%% Calcul via la diagonalisation de la matrice des corrélations
% calcul de la matrice des corrélations
cor = corrcoef(concMat);
% trace de la matrice = inertie totale
trace = sum(diag(cor));
% diagonalisation avec la fonction eig
[V,D,W] = eig(cor);
[d,ind] = sort(diag(D), 'descend');
Ds = D(ind,ind);
eigVect = V(:,ind);
eigVal = sum(Ds); 
Ws = W(:,ind);
% calcul des corrélations des variables avec les composante
corFact = NaN(size(eigVect));
for i = 1:size(concMat,2)
    corFact(:,i) = sqrt(abs(eigVal(i)))*eigVect(:,i);
end

%% Détermination du nombre de composantes à retenir
% On cherche le point d'inflection
% epsilon = diff(eigVal);
% delta = diff(epsilon);
% inflexPt = find(delta < 0, 1);
% 
% figure,
% hold on
% plot(eigVal)
% plot(inflexPt, eigVal(inflexPt), '+')
% hold off
% 
% keptInertia = sum(eigVal(1:inflexPt-1))/sum(eigVal);

% criticalVal = NaN(size(eigVal));
% brStick = 1./(1:length(eigVal));
% for i=1:length(criticalVal)
%     criticalVal(i) = sum(brStick(1:length(eigVal)+1-i));
% end

n = size(concMat,1);
p = size(concMat,2);

%% Analyse des relations entre les variables via les composantes
% figure
% plot(corFact(:,1),corFact(:,2),'.')
% text(corFact(:,1),corFact(:,2),split(num2str(z_)),'FontSize', 5)
% xline(0, '--')
% yline(0, '--')
% circle(0,0,1)

% uiopen(['/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/' ID '-dC.fig'],1)

%% Analyse des proximités entre individus via leurs coordonnées factorielles

[coeff,score,latent] = pca(concMat);

figure(iID), clf,
hold on
% plot(score(:,1),score(:,2))
% plot(score(:,1),score(:,2), '+')

scatter(score(:,1),score(:,2), 10, tC/60/60, 'filled')
xline(0, '--')
yline(0, '--')
hold off
% text(score(:,1),score(:,2),split(num2str(tC/60/60)),'FontSize', 5)
title(ID)
xlim([-2000 1000])
ylim([-200 300])
c = colorbar;
c.Label.String = 'Time (h)';
% lines = get(gca, 'Children');
% set(lines, {'Color'}, nToColorMap(length(lines)))

%     figName = [path runID '-dCstd'];
%     savefig(f3, [figName '.fig']);
%     exportgraphics(f3, [figName '.png']);
    
end, clear iID,

% function h = circle(x,y,r)
% hold on
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% h = plot(xunit, yunit, '--');
% hold off
% end

% branch test