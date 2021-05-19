wind = 1; % km/h
month = 03; 

modSize = 400e-6; % particles size tested (m)
RhoP = 1014; % plage de densite a tester (kg.m⁻³)
nPart = 1e3; % number of particles

tf = 1e5; % simulation time (s)
dt_test = 60*60; % test time interval (s)


%% Data loading
L = 55;
N = 100;
dz= L/N;
z = (0:dz:L)'; % meshes boundaries
z_ = z(1:end-1)+dz/2; % center of the meshes

day = '10fev'; % day corresponding to diffusive turbulence data

NsecTest = 60*60;
[~, ~, ~, ppHist] = varMP_model(modSize, RhoP, false, nPart, tf, dt_test, wind, month, 0, 0,L, N, day,NsecTest);

hConc = NaN(size(ppHist,1),length(z_));
for hStep = 1:size(ppHist,1)
    pp = ppHist(hStep,:);
    if ~isnan(sum(pp))
        histi = histogram(pp, "BinEdges", z, 'Visible', 'off').Values;
        hConc(hStep,:) = histi/dz*L/nPart;
    end
end, clear hStep,
        
meanConc = mean(hConc, 'omitnan');
stdConc = std(hConc, 'omitnan');

figure(1), clf,
plot(meanConc, -z_, meanConc+stdConc, -z_, '--', meanConc-stdConc, -z_, '--')

% séparation des particules en 10 batchs
% nBatch = 1000;
% test = [10 100];
% col = {'r', 'g', 'b', 'm'};
% figure(1), clf, hold on,
% for j=1:length(test)
%     nBatch = test(j);
%     sBatch = nPart/nBatch;
% 
%     iRperm = randperm(length(PartPos));
%     iRperm = reshape(iRperm, [nBatch sBatch]);
% 
%     PartBatch = PartPos(iRperm);
% 
%     % histBatch = zeros(nBatch,length(z_));
%     concBatch = zeros(nBatch,length(z_));
%     for i = 1:nBatch
%         histi = histogram(PartBatch(i,:), "BinEdges", z, 'Visible', 'off').Values;
%         concBatch(i,:) = histi/dz*L/sBatch;
%     end
% 
%     meanConc = mean(concBatch,1)*1e-4; 
%     stdConc = std(concBatch,0,1)*1e-4;
% 
%     p = plot(meanConc, -z_, col{j}, meanConc+stdConc, -z_, ['--' col{j}], meanConc-stdConc, -z_, ['--' col{j}]);
% end
% hold off
