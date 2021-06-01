function [nPartDepth, depthSamples] = getDataNpart(MPtypeName, MPsize, getConcentration, DataFile, SamplingDate)
%%GETDATADONCENTRATION number of particles per sampling depth
% MPtypeName : selected type of particles (string), if all type =false
% MPsize : selected size of particles (m), if all sier =false
% getConcentration : if true, nPartDepth is converted to concentrations
% DataFile : path to the data file
% SamplingDate : date corresponding to the samples of interest 

%% Default arguments
if nargin < 5
    SamplingDate = datetime('3/18/2021');
end
if nargin < 4
    DataFile = '../Data/data_mps.txt';
end
if nargin < 3
    getConcentration = false;
end
if nargin < 2
    MPsize = false;
end
if nargin == 0
    MPtypeName = false;
end
    
%% Load data file
dataTable = load_MPs_data(DataFile);


%% Filter Data
% Sampling day
condition = dataTable(:,'date').Variables == SamplingDate; 
% Type and size
if MPtypeName
    typeSTR =  containers.Map({'fibre','fragment','film','mousse'},0:3); % List of existing type
    MPtype = typeSTR(MPtypeName);
    condition = condition & dataTable(:,'type').Variables == MPtype;
end
if MPsize
    condition = condition & dataTable(:,'med_size').Variables == MPsize*1e6;
end
    
filteredDataTable = dataTable(condition,:); 

%% Sum the number of particles at each depth

depthSamples = unique(filteredDataTable(:,'depth').Variables); % sorted list of depth

% Find number of particles at each depth
nPartDepth = zeros(size(depthSamples)); % preallocation 
for i=1:length(depthSamples)
    cond = filteredDataTable(:,'depth').Variables == depthSamples(i);
    nPartDepth(i) = sum(filteredDataTable(cond,'n').Variables);
end, clear i,

if getConcentration
    % Compute concentration at each depth
    VolumeList = [88.78314375;331;435;437;392;291]; % Volume list for default data/date NEEDS TO BE CHANGED
    a = nPartDepth;
    nPartDepth = a./VolumeList*8;
end

end 