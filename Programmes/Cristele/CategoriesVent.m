ModeleHydro='2012RHOMA_arome_003.nc';
indTf=60*24;indTd=1;
AfficheDessin=false;
AfficheDessin=true;
%ncdisp('2012RHOMA_arome_003.nc');

%nc=netcdf(ModeleHydro);
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)
load Diffusion.mat
Lon0= 5.29;Lat0=43.24; %point Somlit

% [I0,J0,D0]=ReperePoint(Lon,Lat,Lon0,Lat0);

% TX=squeeze(nc{'TAUX'}(:,:,I0,J0));
% TY=squeeze(nc{'TAUY'}(:,:,I0,J0));
% KZ_=squeeze(nc{'KZ'}(:,:,I0,J0));
% close(nc)

TAU=sqrt(TauX.^2+TauY.^2);
TAU_=TAU(indTd:indTf);
TAU_=TAU;
Temps_=Temps(indTd:indTf);

%Force 0 [0 ; 1 km/h[ 
%Force 1 [1 ; 5]      
%Force 2 [6 ; 11]
%Force 3 [12 ; 19]
%Force 4 [20 ; 28]
%Force 5 [29 ; 38]
%Force 6 [39 ; 49]
%Force 7 [50 ; 61]
%Force 8 [62 ; 74]
%Force 9 [75 ; 88]
%Force 10 [89 ; 102]
%Force 11 [103 ; 117]
%Force 12 [118 ; +inf[

SeuilsVitessesKMH= [0 1 6 12 20 29 39 50 62 75 89 103 118]; %km/h
%SeuilsVitessesKMH= [7.5:11.5]; %km/h
SeuilsVitessesMS= SeuilsVitessesKMH/3.6; %m/s

%  CdSeuils=10^-3*(0.63+0.066*SeuilsVitessesMS); %Smith & Banke
CdSeuils=10^-3*(0.43+0.096*SeuilsVitessesMS); %Geernaert

    
SeuilsTAU=1.292*CdSeuils.*SeuilsVitessesMS.*SeuilsVitessesMS;
% plot(SeuilsTAU,SeuilsVitessesKMH,'.')
Force=NaN*TAU_;
for indice=1:(size(SeuilsTAU,2)-1)
    ii=find(TAU_<SeuilsTAU(indice+1) & TAU_>=SeuilsTAU(indice));
    if size(ii,1)~=0
        Force(ii)=indice-1;
        KZmoy(indice,:)=mean(KZ0(ii,:),1);
        KZmed(indice,:)=median(KZ0(ii,:),1);
%         KZmed_(indice,:)=KZ0(ii(floor(size(ii,1)/2),:);        
        if AfficheDessin
            figure,clf
            plot(KZ0(ii,:),z0,'k',KZmoy(indice,:),z0,'r',KZmed(indice,:),z0,'g')
            title(['Force : ' num2str(indice-1)])
        end
    end
end
% if AfficheDessin
%     figure(2),clf
%     plot(Temps_-T0,TAU_),hold on
%     scatter(Temps_-T0,TAU_,50,Force,'filled')
% end


% figure(2)
% plot(KZmed,z0,'--')
% legend('Force0','Force1','Force2','Force3','Force4','Force5','Force6','Force7','Force8','Force9','Location', 'best')
% title('Categories de Vent Median')
% 
% hold on
% 
% plot(KZmoy,z0,'-')
% legend('Force0','Force1','Force2','Force3','Force4','Force5','Force6','Force7','Force8','Force9','Location', 'best')
% title('Categories de Vent Moyen')

