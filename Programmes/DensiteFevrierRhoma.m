AfficheDessin=false;
GainTemps=1;
FichDiffusion='/media/ian/Elements/Ian_Plastique/Data/Diffusion';
FichSomlit='/media/ian/Elements/Ian_Plastique/Data/DataSomlit';
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];

load(SauvegardeModeleHydro)
load(FichDiffusion)
TAU=sqrt(TauX.^2 + TauY.^2);

% plot(TAU)

TAUJanFev=TAU(1:24*(31+29));

%%%Contampump 10.02.20 (10:22 a 11:54 Vent moyen=7.7 Force2)
Vitesse10FevKH=7.7;
Vitesse10FevMS=Vitesse10FevKH/3.6;

Cd10Fev=10^-3*(0.43+0.096*Vitesse10FevMS); %Geernaert

TAU10Fev=1.292*Cd10Fev.*Vitesse10FevMS.*Vitesse10FevMS;

% RHOMAJTAU10=find(TAU==TAU10Fev);
JRhomaTAU10=find(abs(TAUJanFev-TAU10Fev) < 0.0001);

%%%Contampump 03.02.20 (10:30 a 11:11 Vent moyen=10.05 Force2)
Vitesse03FevKH=10.05;
Vitesse03FevMS=Vitesse03FevKH/3.6;

Cd03Fev=10^-3*(0.43+0.096*Vitesse03FevMS); %Geernaert

TAU03Fev=1.292*Cd03Fev.*Vitesse03FevMS.*Vitesse03FevMS;

JRhomaTAU03=find(abs(TAUJanFev-TAU03Fev) < 0.000005);


%%%--------------Donnees Temperature/Salinite Rhoma 10 Fevrier

Lon0= 5.29;Lat0=43.24; %point Somlit
Temps0=JRhomaTAU10/24; 
% ncdisp('2012RHOMA_arome_003.nc')
T0=datenum(str2num(datestr(Temps(1),'yyyy')),1,1)-1;

[tt,iT0]=min(abs(Temps-T0-Temps0));

if AfficheDessin
    figure(1),clf,hold on
    plot(Temps-T0,sqrt(TauX.^2+TauY.^2))
    plot(Temps0,sqrt(TauX(iT0).^2+TauY(iT0).^2),'*r')
    xlabel('Temps - jour Julien')
    ylabel('|Stress du au vent|')
    title(['jour ',num2str(Temps0),' - Point SOMLIT'])
end

%[I0,J0,D0]=ReperePoint(Lon,Lat,Lon0,Lat0);

    KZ_Fev10=KZ0(iT0,:);
    Sal_Fev10=Salinite0(iT0,:);
    Temp_Fev10=Temperature0(iT0,:);
    z_Fev10=z0;
    z__Fev10=z0_;

    Row_Fev10=CalculDensite(Temp_Fev10,Sal_Fev10); %source: edu.obs-mip.fr

if AfficheDessin
    figure(2),clf
    minz0_=min(z__Fev10);
    minz0=min(z__Fev10);
    maxT=max(Temp_Fev10);
    minT=min(Temp_Fev10);
    maxsal=max(Sal_Fev10);
    minsal=min(Sal_Fev10);
    minKZ0=min(KZ_Fev03);
    maxKZ0=max(KZ_Fev03);

    subplot(1,3,1),plot(Temp_Fev10,z__Fev10),xlabel('T'),ylabel('Depth (m)')
    axis([minT maxT minz0_ 0])
    subplot(1,3,2),plot(Sal_Fev10,z__Fev10),xlabel('Sal')
    axis([minsal maxsal minz0_ 0])
    title(['jour ',num2str(Temps0),' - Point SOMLIT'])
    subplot(1,3,3),plot(KZ_Fev10,z_Fev10),xlabel('Ks')
    axis([minKZ0 maxKZ0 minz0 0])

    figure(3),clf
    maxrow=max(Row_Fev10);
    minrow=min(Row_Fev10);
    plot(Row_Fev10,z__Fev10),xlabel('Seawater density (kg.m^-^3)'), ylabel('Depth (m)')
    title(['jour ',num2str(Temps0),' - Point SOMLIT'])
    axis([minrow maxrow minz0 0])
end

%%%------------Donnees Temperature/Salinite Rhoma 03 Fevrier

Lon0= 5.29;Lat0=43.24; %point Somlit
Temps0=JRhomaTAU03/24; 
ModeleHydro='2012RHOMA_arome_003.nc'
% ncdisp('2012RHOMA_arome_003.nc')
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)
T0=datenum(str2num(datestr(Temps(1),'yyyy')),1,1)-1;

load(SauvegardeModeleHydro)

[tt,iT0]=min(abs(Temps-T0-Temps0));

if AfficheDessin
    figure(1),clf,hold on
    plot(Temps-T0,sqrt(TauX.^2+TauY.^2))
    plot(Temps0,sqrt(TauX(iT0).^2+TauY(iT0).^2),'*r')
    xlabel('Temps - jour Julien')
    ylabel('|Stress du au vent|')
    title(['jour ',num2str(Temps0),' - Point SOMLIT'])
end

%[I0,J0,D0]=ReperePoint(Lon,Lat,Lon0,Lat0);

    KZ_Fev03=KZ0(iT0,:);
    Sal_Fev03=Salinite0(iT0,:);
    Temp_Fev03=Temperature0(iT0,:);
    z_Fev10=z0;
    z__Fev10=z0_;

    Row_Fev03=CalculDensite(Temp_Fev03,Sal_Fev03); %source: edu.obs-mip.fr

if AfficheDessin
    figure(4),clf
    minz0_=min(z__Fev03);
    minz0=min(z__Fev03);
    maxT=max(Temp_Fev03);
    minT=min(Temp_Fev03);
    maxsal=max(Sal_Fev03);
    minsal=min(Sal_Fev03);
    minKZ0=min(KZ_Fev03);
    maxKZ0=max(KZ_Fev03);

    subplot(1,3,1),plot(Temp_Fev03,z__Fev03),xlabel('T'),ylabel('Depth (m)')
    axis([minT maxT minz0_ 0])
    subplot(1,3,2),plot(Sal_Fev03,z__Fev03),xlabel('Sal')
    axis([minsal maxsal minz0_ 0])
    title(['jour ',num2str(Temps0),' - Point SOMLIT'])
    subplot(1,3,3),plot(KZ_Fev03,z_Fev03),xlabel('Ks')
    axis([minKZ0 maxKZ0 minz0 0])

    figure(5),clf
    maxrow=max(Row_Fev03);
    minrow=min(Row_Fev03);
    plot(Row_Fev03,z__Fev03),xlabel('Seawater density (kg.m^-^3)'), ylabel('Depth (m)')
    title(['jour ',num2str(Temps0),' - Point SOMLIT'])
    axis([minrow maxrow minz0 0])
end


