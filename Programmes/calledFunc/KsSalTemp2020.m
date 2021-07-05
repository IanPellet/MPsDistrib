function [KZ_day,Row_day,z_day,z__day,Sal_day,Temp_day] = KsSalTemp2020(WindSpeed_kmh, date)
%KSSALTEMP Summary of this function goes here
% Find the day in the 2012RHOMA_arome database with the closest wind at RN2
% situation to return the values of diffusivity, salinity, temperature and
% water density that day.
%   WindSpeed_kmh : Wind Speed in km/h
%   date : (datetime) sampling date

% clear
% date = datetime('18/03/2012','InputFormat', 'dd/MM/yyyy');
% WindSpeed_kmh = 50;
date = datetime(2020,month(date),day(date));
monthStart = date - calmonths(2);
monthEnd = date + calmonths(2);
origin = datetime(0,1,1,0,0,0);
tStart = days(days(monthStart - origin));
tEnd = days(days(monthEnd - origin));
tDate = days(days(date - origin));

load("waterCol_RN2_2020.mat", 'Sal', 'Temp', 'Kz', 'z0', 'TauX', 'TauY', 't');
t = days(t);

TAU=sqrt(TauX.^2 + TauY.^2); % TAU = tension de surface liée au vent

% Extract tau on the 2 month period around month inputed
itStart = find(abs(t-tStart) == min(abs(t-tStart)));
itEnd = find(abs(t-tEnd) == min(abs(t-tEnd)));
itDate = find(abs(t-tDate) == min(abs(t-tDate)));
TAU_period = TAU(itStart:itEnd); % Extract data from dStart to dStop

%%%Contampump 10.02.20 (10:22 a 11:54 Vent moyen=7.7 Force2)
WindSpeed_ms = WindSpeed_kmh/3.6; % Wind speed m/s

Cd_day=10^-3*(0.43+0.096*WindSpeed_ms); %Geernaert               %% ????

TAU_day=1.292*Cd_day.*WindSpeed_ms.*WindSpeed_ms;

% Find a day with similar wind in our data
iTAU_period =find(abs(TAU_period-TAU_day) == min(abs(TAU_period-TAU_day))); % index

% If several days are found, we take the one closer to the middle of the
% month
if length(iTAU_period ) > 1   
    dDist = abs(iTAU_period-itDate);
    dClose = min(dDist);
    iTAU_period = iTAU_period(dDist == dClose);
end

iDay = find(TAU == TAU_period(iTAU_period)); % day index


KZ_day = Kz(iDay,:);
Sal_day = Sal(iDay,:);
Temp_day = Temp(iDay,:);
z_day = z0(iDay,:);
z__day = z_day(1:end-1) + diff(z_day(:))'/2;
 
Row_day = CalculDensite(Temp_day,Sal_day); %source: edu.obs-mip.fr


end

