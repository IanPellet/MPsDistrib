function [CPartDepth, depthSamples] = getDataConcentration(DataFile, VolumeList, MPtype, MPsize)

% Default arguments
if nargin == 0
    DataFile = '../Data/data_mps.txt';
%     VolumeData = '../Data/volumesFiltres.txt';
    VolumeList = [88.78314375;331;435;437;392;291];
    SamplingDate = datetime('3/18/2021');
    MPtype = 0;
    MPsize = 750;
end

dataTable = load_MPs_data(DataFile); % load data file
% volumeTable = readtable(VolumeData); % load volumes file
conditionDate = dataTable(:,'date').Variables == SamplingDate;
depthSamples = unique(dataTable(conditionDate,'depth').Variables); % sorted list of depth

% Filter data for type and size
condition = dataTable(:,'type').Variables == MPtype &...
    dataTable(:,'med_size').Variables == MPsize;
    
filteredDataTable = dataTable(condition,:); 

% Find number of particles at each depth
nPartDepth = zeros(size(depthSamples)); % preallocation 
for i=1:length(depthSamples)
    cond = filteredDataTable(:,'depth').Variables == depthSamples(i);
    nPartDepth(i) = sum(filteredDataTable(cond,'n').Variables);
end, clear i,

% Compute concentration at each depth
CPartDepth = nPartDepth./VolumeList;

end 