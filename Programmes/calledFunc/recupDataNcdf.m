function [fileName] = recupDataNcdf(station)
%% RECUPDATANCDF Saves data at station's point 
% Saves salinity, temperature, diffusivity, depth, data points levels
% at the data point corresponding to the station and the time array t.
%
% station : (string) station's name
%
% Data saved as '../Data/waterCol_station.mat'

start % add necessary paths 

% Load file with indexes corresponding to each stations
stationFile = '../Data/stationIJ_CEREGE.mat';
load(stationFile,'stationIJ');
% Get the indexes of the right station
I0 = stationIJ{stationIJ{:,'station'} == station,'I0'};
J0 = stationIJ{stationIJ{:,'station'} == station,'J0'};

% Load netCDF data files
ncSal = netcdf('../../rhoma2012/SAL.nc');
ncTemp = netcdf('../../rhoma2012/TEMP.nc');
ncKz = netcdf('../../rhoma2012/KZ.nc');

% Get the full data arrays pointers
fullSal = ncSal{'SAL'};
fullTemp = ncTemp{'TEMP'};
fullKz = ncKz{'KZ'};

% Get the data at the station's point
Sal = fullSal(:,:,I0,J0);
Temp = fullTemp(:,:,I0,J0);
Kz = fullKz(:,:,I0,J0);
% Get time array
t = seconds(ncKz{'time'}(:,:));

H0 = ncKz{'H0'}(I0,J0); % bathymetry relative to the mean level
eta = ncKz{'XE'}(:,I0,J0); % sea surface height
sigma_w = ncKz{'level_w'}(:,:); % sigma level at the interface

% Compute depth of data points
z0 = zeros(length(t), length(sigma_w));
for i=1:length(t)
    z0(i,:) = eta(i) + sigma_w*(H0+eta(i));
end

% Save data to file
fileName = join(['../Data/waterCol_' station ".mat"],"");
save(fileName, 'Sal', 'Temp', 'Kz', 't', 'H0', 'z0');

disp(join(['Saved to' fileName]))
end
