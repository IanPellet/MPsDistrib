path = '../Results/varMP/';

wind = 10; % km/h
month = 02; % February

load('../Data/stationLonLat_CEREGE.mat','station') % load position data of stations
staName = 'RN1'; % chosen station
iStation = find(station(:,'Station').Variables == staName); % corresponding row of the table station
% Position of station
lon0 = station(1,'Lon').Variables;
lat0 = station(1,'Lat').Variables;

depth_res = [0 1]; % depth of samples (m)

size_test = [1250 1750]*1e-6; % particles size tested (m)
% rhop_test = linspace(850,1300, 10);
rhop_test = [1010.5]; % particles density (kg.m⁻³)
type_dict = containers.Map({'fibre','fragment','film','mousse'},0:3); % possible types
type_name = 'fragment'; % choose a type
type = type_dict(type_name); % corresponding int

nPart = 50e3; % number of particles

tf = 1e5; % simulation time (s)
dt_test = 60*60; % test time interval (s)

% Resultats = cell(length(D_test)*length(rhop_test),4); % Memory allocation for storing results
clear Resultats

i = 0;
for D = size_test
    for rhop = rhop_test
        i = i+1;
        [Zi, Ci] = varMP_model(D, rhop, type, nPart, tf, dt_test, wind, month, lon0, lat0, path);
%         Resultats(i,:) = {D, rhop, Zi, Ci}; % Store values
        Resultats(i).D = D;
        Resultats(i).rhop = rhop;
        Resultats(i).z = Zi;
        Resultats(i).C = Ci;
        Resultats(i).depth = depth_res;
        
        Ctest = zeros(depth_res);
        for j = 1:length(depth_res)
            depth = depth_res(j);
            deltat_depth = abs(Zi-depth);
            deltat_min = min(deltat_depth);
            id = find(deltat_depth == deltat_min,1,'first');
            Cd = Ci(id);
            
            Ctest(j) = Cd;
        end
        
        Resultats(i).C_depth = Ctest;
    end
end

file_name = ['results_D',num2str(min(D)),'-',num2str(length(D)),'-',num2str(max(D)),...
    '_rhop',num2str(rhop), '_nPart',num2str(nPart), '_tf', num2str(tf),...
    '_dtest', num2str(dt_test),'_wind', num2str(wind), '_month', num2str(month)];
save([path,file_name,'.mat'],'Resultats');


