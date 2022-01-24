function [allConcMat, allTauxMat, MPConcMat, AggConcMat] = plotAggPcolor(IDrun, N, dz, nMP, nAgg)
filePath = '/media/ian/Transcend/MPsDistrib/Results/AggPlot/';

nbrFiles = GetiFiles(IDrun);


allConc = cell(length(nbrFiles),1);
MPConc = cell(length(nbrFiles),1);
AggConc = cell(length(nbrFiles),1);
allTaux = cell(length(nbrFiles),1);
for iN = 1:length(nbrFiles)
    load([filePath IDrun '_' num2str(nbrFiles(iN)) '.mat'])
    
    while isempty(zHistory{end})
        zHistory = zHistory(1:end-1);
    end
    
    allConc{iN} = GetConc(zHistory, 1, nMP+nAgg);
    MPConc{iN} = GetConc(zHistory, 1, nMP);
    AggConc{iN} = GetConc(zHistory, nMP+1, nMP+nAgg);
    
    while isnan(tauxAgg(end))
        tauxAgg = tauxAgg(1:end-1);
    end
    allTaux{iN} = tauxAgg;
end, clear iN zHistory tauxAgg hConc,
allConcMat = cell2mat(allConc);
MPConcMat = cell2mat(MPConc);
AggConcMat = cell2mat(AggConc);
allTauxMat = cell2mat(allTaux);

save([filePath 'allConc_' IDrun '.mat'], 'allConcMat', 'allTauxMat', 'MPConcMat', 'AggConcMat')
% 
% %% Plot
% h = pcolor(allConcMat');
% cb = colorbar;
% % caxis([0 max(max(meanC))])
% set(h, 'EdgeColor', 'none');
end