function [Keep] = KeepSizeType_threshold(DataFile, SamplingDate, threshold)
%% Find size and type of particles above a number of particle defined by threshold 

if nargin == 0
    DataFile = '../Data/data_mps.txt'; % File containing data
    SamplingDate = datetime('3/18/2021'); % Date of the data of interest
    threshold = 10; % Minimal number of particles
end

dataTable = load_MPs_data(DataFile); % load data file
% Only keep the data of the given date
conditionDate = dataTable(:,'date').Variables == SamplingDate;
dataMultinet = dataTable(conditionDate,:);

sizeSamples = unique(dataMultinet(:,'med_size').Variables); % sorted list of particle size
typeSamples = unique(dataMultinet(:,'type').Variables); % sorted list of particle type


NumbPart = zeros(length(typeSamples),length(sizeSamples)); % memory allocation for number of particles of a given type and size
for i=1:length(typeSamples)
    for j = 1:length(sizeSamples)
        MPsize = sizeSamples(j);
        MPtype = typeSamples(i);
        
        % Filter data for type and size
        condition = dataMultinet(:,'type').Variables == MPtype &...
        dataMultinet(:,'med_size').Variables == MPsize;
        filteredDataTable = dataMultinet(condition,:); 

        NumbPart(i,j) = sum(filteredDataTable(:,'n').Variables); % number of particles of type MPtype and size MPsize
    end, clear j,
end, clear i MPsize MPtype condition filteredDataTable,

typeSTR = {'fibre','fragment','film','mousse'}; % List of existing type

% Plot number of particles of each size for each type
figure(1), clf,
plot(sizeSamples,ones(size(sizeSamples))*threshold,'--', 'DisplayName', 'Threshold')
hold on
for i=1:length(typeSamples)
    plot(sizeSamples,NumbPart(i,:),'h','DisplayName',typeSTR{typeSamples(i)+1})
end, clear i,
xticks(sizeSamples)
xlabel('Particle size (m)')
ylabel('Number of particles')
legend('Location', 'best')
hold off

% Find the row and column indexes of data above the threshold
[row,col] = find(NumbPart >= threshold);
Keep = struct('Type', typeSTR(row), 'Size', num2cell(sizeSamples(col))');

clear conditionDate DataFile dataMultinet dataTable depthSamples filteredDataTable SamplingDate sizeSamples typeSamples typeSTR threshold row col keep
end