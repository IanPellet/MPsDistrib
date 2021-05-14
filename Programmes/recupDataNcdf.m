clear

station = "RN2";
date = datetime('18/03/2012','InputFormat', 'dd/MM/yyyy');
monthStart = date - calmonths(1);
monthEnd = date + calmonths(1);
origin = datetime(1900,1,1,0,0,0);
tStart = monthStart - origin;
tEnd = monthEnd - origin;

stationFile = '../Data/stationIJ_CEREGE.mat';
load(stationFile,'stationIJ');
I0 = stationIJ{stationIJ{:,'station'} == station,'I0'};
J0 = stationIJ{stationIJ{:,'station'} == station,'J0'};

ncSal = netcdf('../../rhoma2012/SAL.nc');
ncTemp = netcdf('../../rhoma2012/TEMP.nc');
ncKz = netcdf('../../rhoma2012/KZ.nc');

Sal = ncSal{'SAL'};
Temp = ncTemp{'TEMP'};
Kz = ncKz{'KZ'};
t = seconds(ncKz{'time'}(:,:));

itStart = find(abs(t-tStart) == min(abs(t-tStart)));
itEnd = find(abs(t-tEnd) == min(abs(t-tStart)));


Kz_month = Kz(itStart:itEnd,:,I0,J0);
