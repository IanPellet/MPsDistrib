% netCDF file handling


clear,

sta = "RN2";
% Load file with lon/lat corresponding to each stations
stationFile = '../Data/stationLonLat_CEREGE.mat';
load(stationFile,'station');
% Get the lon/lat of the right station
Lon = station{strcmp(station{:,'Station'}, sta),'Lon'};
Lat = station{strcmp(station{:,'Station'}, sta),'Lat'};

clear sta stationFile station,

% Load netCDF data file
ncFile = netcdf('../Data/rhoma2020/2020RHOMA_WRF6h_003.nc');

longitude_u = ncFile{'longitude_u'}(:,:); % longitude at u location
longitude_v = ncFile{'longitude_v'}(:,:); % longitude at v location
latitude_u = ncFile{'latitude_u'}(:,:); % latitude at u location
latitude_v = ncFile{'latitude_v'}(:,:); % latitude at v location

[~,i_u] = min(abs(unique(longitude_u)-Lon));
[~,i_v] = min(abs(unique(longitude_v)-Lon));
[~,j_u] = min(abs(unique(latitude_u)-Lat));
[~,j_v] = min(abs(unique(latitude_v)-Lat));

clear longitude_u longitude_v latitude_u latitude_v Lon Lat 

t = seconds(ncFile{'time'}(:,:));
date = datetime(2020,2,13,7,0,0);
origin = datetime(1900,1,1,0,0,0);
tDate = seconds(seconds(date - origin));
iDay = find(abs(t-tDate) == min(abs(t-tDate)),1);

clear t date origin tDate,

UZ = ncFile{'UZ'}(iDay,:,:,:); % sea_water_x_velocity_at_u_location
VZ = ncFile{'VZ'}(iDay,:,:,:); % sea_water_y_velocity_at_v_location

dx_u = ncFile{'dx_u'}(:,:); % mesh size along x at u location
dx_v = ncFile{'dx_v'}(:,:); % mesh size along x at v location
dy_u = ncFile{'dy_u'}(:,:); % mesh size along y at u location
dy_v = ncFile{'dy_v'}(:,:); % mesh size along y at v location

dudx = deriveSpeed(UZ,dx_u,1);
dudy = deriveSpeed(UZ,dy_u,2);
dvdx = deriveSpeed(VZ,dx_v,1);
dvdy = deriveSpeed(VZ,dy_v,2);

clear UZ VZ dx_u dx_v dy_u dy_v,




% 
% 
% % Get the data at the station's point
% % ni_u = ncFile{'ni_u'}(iDay,:,:,:); % x-dimension of the grid at u location
% % nj_u = ncFile{'nj_u'}(iDay,:,:,:); % y-dimension of the grid at u location
% % ni_v = ncFile{'ni_v'}(iDay,:,:,:); % x-dimension of the grid at v location
% % nj_v = ncFile{'nj_v'}(iDay,:,:,:); % y-dimension of the grid at v location
% 
% 
% SIG = ncFile{'SIG'}(:,I0,J0); % sigma variable
% 
% 
% % Get time array
% t = seconds(ncFile{'time'}(:,:));

function dudx = deriveSpeed(UZ,dx_u,dim)
    dx_u_rs = repmat(reshape(dx_u,[1 size(dx_u)]),[size(UZ,1),1,1]);
    if dim == 1
        dx_u_rs_ = (dx_u_rs(:,1:end-1,:)+dx_u_rs(:,2:end,:))/2;
    elseif dim == 2
        dx_u_rs_ = (dx_u_rs(:,:,1:end-1)+dx_u_rs(:,:,2:end))/2;
    end
    dudx = diff(UZ,1,dim+1)./dx_u_rs_;
end