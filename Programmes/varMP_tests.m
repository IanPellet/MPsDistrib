path = '../Results/varMP/';

load('../Data/stationLonLat_CEREGE.mat','station')

wind = 10; % km/h
month = 02;
lon0 = station(1,'Lon').Variables;
lat0 = station(1,'Lat').Variables;

depth_res = [0 1]; % m

D_test = [1250 1750]*1e-6;
% rhop_test = linspace(850,1300, 10);
rhop_test = [1010.5]; % kg.m⁻³
nPart = 50e3;

dt = 10;
tf = 1e5;
dt_test = 60*60;

% Resultats = cell(length(D_test)*length(rhop_test),4); % Memory allocation for storing results
clear Resultats

i = 0;
for D = D_test
    for rhop = rhop_test
        i = i+1;
        [Zi, Ci] = varMP_model(D, rhop, nPart, dt, tf, dt_test, wind, month, lon0, lat0, path);
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
            id = find(deltat_depth == deltat_min);
            Cd = Ci(id);
            
            Ctest(j) = Cd;
        end
        
        Resultats(i).C_depth = Ctest;
    end
end

file_name = ['results_D',num2str(min(D)),'-',num2str(length(D)),'-',num2str(max(D)), '_rhop',num2str(rhop),...
    '_nPart',num2str(nPart), '_dt', num2str(dt),...
    '_tf', num2str(tf), '_dtest', num2str(dt_test),'_wind', num2str(wind),...
    '_month', num2str(month)];
save([path,file_name,'.mat'],'Resultats');


