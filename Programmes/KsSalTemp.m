function [KZ_day,Row_day,z_day,z__day] = KsSalTemp(WindSpeed_kmh, month)
%KSSALTEMP Summary of this function goes here
% Find the day in the 2012RHOMA_arome database with the closest wind
% situation to return the values of turbidity, salinity, temperature and
% water density that day.
%   WindSpeed_kmh Wind Speed in km/h

ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)
% Loaded variables :  Adv - C - C_ - Diff - H0 - ii - Lat - Lon - Nu -
% Sigma - TauX - TauY - Temps 

FichDiffusion='/media/ian/Elements/Ian_Plastique/Data/Diffusion';
load(FichDiffusion)
% Loaded variables :  KZ0 - Row0 - Salinite0 - Temperature0 - z0 - z0_

mDays = [31 29 31 30 31 30 31 31 30 31 30 31]; % Days in each month

TAU=sqrt(TauX.^2 + TauY.^2); % TAU = tension de surface liée au vent

%%% NE GERE PAS LES CAS SPECIAUX DES MOIS DE JANVIER ET DECEMBRE %%%

% Extract tau on the 3 month period around month inputed
dStart = 24*sum(mDays(1:month-1)); % start one month before
dStop = 24*sum(mDays(1:month+1)); % stop one month after
TAU_period = TAU(dStart:dStop); % Extract data from dStart to dStop

%%%Contampump 10.02.20 (10:22 a 11:54 Vent moyen=7.7 Force2)
WindSpeed_ms = WindSpeed_kmh/3.6; % Wind speed m/s

Cd10Fev=10^-3*(0.43+0.096*WindSpeed_ms); %Geernaert               %% ????

TAU10Fev=1.292*Cd10Fev.*WindSpeed_ms.*WindSpeed_ms;

JRhomaTAU =find(abs(TAU_period-TAU10Fev) < 0.0001); % Find a day with similar wind in our data

% If several days are found, we take the one closer to the middle of the
% month
if length(JRhomaTAU ) > 1
    dMid = mean(dStart,dStop);
    dDist = abs(JRhomaTAU-dMid);
    dClose = min(dDist);
    JRhomaTAU = JRhomaTAU(dDist == dClose);
end

%%%--------------Donnees Temperature/Salinite Rhoma 10 Fevrier

Temps0=JRhomaTAU/24; 
T0=datenum(str2num(datestr(Temps(1),'yyyy')),1,1)-1;

[tt,iT0]=min(abs(Temps-T0-Temps0));


KZ_day = KZ0(iT0,:);
Sal_day = Salinite0(iT0,:);
Temp_day = Temperature0(iT0,:);
z_day = z0;
z__day = z0_;
 
Row_day = CalculDensite(Temp_day,Sal_day); %source: edu.obs-mip.fr




end

