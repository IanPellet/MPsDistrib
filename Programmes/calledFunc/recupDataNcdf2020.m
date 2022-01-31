function [fileName] = recupDataNcdf2020(station)
%% RECUPDATANCDF2020 Saves data at station's point 
% Saves salinity, temperature, diffusivity, depth, data points levels
% at the data point corresponding to the station and the time array t.
%
% station : (string) station's name
%
% Data saved as '../Data/waterCol_station.mat'

% start % add necessary paths 

% ncdisp("../Data/rhoma2020/2020RHOMA_WRF6h_003.nc",'/','min')


% Load file with indexes corresponding to each stations
stationFile = '../Data/stationIJ_CEREGE.mat';
load(stationFile,'stationIJ');
% Get the indexes of the right station
I0 = stationIJ{stationIJ{:,'station'} == station,'I0'};
J0 = stationIJ{stationIJ{:,'station'} == station,'J0'};

% Load netCDF data file
ncFile = netcdf('../Data/rhoma2020/2020RHOMA_WRF6h_003.nc');

% Get the full data arrays pointers
% fullSal = ncFile{'SAL'};
% fullTemp = ncFile{'TEMP'};
% fullKz = ncFile{'KZ'};

<<<<<<< HEAD

% Get the data at the station's point
Sal = ncFile{'SAL'}(:,:,I0,J0);
Temp = ncFile{'TEMP'}(:,:,I0,J0);
Kz = ncFile{'KZ'}(:,:,I0,J0);
TauX = ncFile{'TAUX'}(:,I0,J0);
TauY = ncFile{'TAUY'}(:,I0,J0);

% Get time array
t = seconds(ncFile{'time'}(:,:));


H0 = ncFile{'H0'}(I0,J0); % bathymetry relative to the mean level
eta = ncFile{'XE'}(:,I0,J0); % sea surface height
sigma_w = ncFile{'level_w'}(:,:); % sigma level at the interface

% Uz = ncFile{'UZ'}(:,:,I0,J0);
% Vz = ncFile{'VZ'}(:,:,I0,J0);
% Chl = ncFile{'phy_phyto_Chl'}(:,:,I0,J0);

% Compute depth of data points
z0 = zeros(length(t), length(sigma_w));
for i=1:length(t)
    z0(i,:) = eta(i) + sigma_w*(H0+eta(i));
end
toc
% Save data to file
fileName = join(['../Data/waterCol_' station "_2020.mat"],"");
% save(fileName, 'Sal', 'Temp', 'Kz', 't', 'H0', 'z0', 'TauX', 'TauY', 'Uz', 'Vz', 'Chl');

save(fileName, 'Sal', 'Temp', 'Kz', 't', 'H0', 'z0', 'TauX', 'TauY');

disp(join(['Saved to' fileName]))
end
