% Ws qui minimise l'Err ...
% (Csurf peut etre extrapole par C(-1m) )
% 
% On connaît :
%     dx
%     C_
%     Z_
%     dh
global dt dx

Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];

% Initilisation :
D=350e-6; %m : Diametre
Lon0= 5.29;Lat0=43.24; %point Somlit
ModeleHydro='2012RHOMA_arome_003.nc';
indNom=3;

SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)
%[I0,J0]=ReperePoint(Lon,Lat,Lon0,Lat0);

tf= 100*86400; dtmax=0.01; 
Tdes=60*60; dConcMax=5E-5; % Tdes : intervalle de temps entre les tests d'équilibre 
% dConcMax : seuil de delta de concentration à partir duquel on de considère à l'équilibre 
dh=0.15; % profondeur sur laquelle le filet prélève

% Eulerien vs Lagrangien
tf= 100*86400; % maximum simulation time
dt_max=0.01; % maximun time interval
dt_test = 60*30; % time interval between equilirium tests

clear Concentration err
%L=H0(I0,J0);   %m
L = 50; % Profondeur (arbitraire)
N=2000;  dx= L/N;  z=0:dx:L; % x : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % milieu de chaque maille
Xmin=0;Xmax=L;Cmin=-1;Cmax=2;
        %z0=H0(I0,J0)*Sigma;
        z0=L*Sigma;
        %z0_=H0(I0,J0)*(Sigma(1:end-1)+Sigma(2:end))/2;
        z0_=L*(Sigma(1:end-1)+Sigma(2:end))/2;    

% conditions initiales 
CMes=[0.62 0.34 0.06 0.02 0]; % concentrations mesurées
ZMes=[1 10 15 40 L]; % profondeur de chaque mesure
C = interp1(ZMes,CMes,z(1:end-1)+dx/2,'pchip'); % interpolation sur x
C=max(0*C,C); C_init = C;

figure(1),clf,plot(C,-z_,'r',CMes,-ZMes,'og'),xlabel('Concentration (kg.m^-^3)'), ylabel('Depth (m)')


    
%% 1) calculer le profil de concentration associé à Ws
Concentration(1,:)=C;
row = 1000;

% Determiner rho eau
DensiteFevrierRhoma
%Nu=interp1(z0,KZ_Fev10,-x_,'pchip');
Ks = 0.01;
Nu = ones(1,2000)*Ks;

figure(4),
    subplot(1,2,1), plot(row,-z)
    subplot(1,2,2), plot(Nu,-z_)

InitialisationVitesseTransport

rop=1010.5;
S=rop./row;     D_=((g*(abs(S-1))/nuw^2).^(1/3))*D;
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
u=Ws; u(rop<row)=-Ws(rop<row);

% Analytical solution at equilibrium
%Ccalc = C_analytical(Ws, Ks,x_, C);

u0_=max(u);Nu0_=max(Nu);
if u0_~=0 & Nu0_~=0; 
   dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
elseif u0_==0 & Nu0_~=0; 
   dt=dx*dx/(2*Nu0_)*0.5;
elseif u0_~=0 & Nu0_==0; 
   dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
else, 
   dt=dtmax;
end

t=0; OnContinue=true;
while OnContinue
   t=t+dt;
   C=StepTransport (u,Nu,C,'UpWind'); 

   if (mod(t,Tdes)<=dt/2 | Tdes-mod(t,Tdes)<=dt/2 )
       index=round(t/Tdes)+1;

       Concentration(index,:)=C;
       Ecart(index,:)=max((Concentration(index-1)-Concentration(index)).^2);
       if Ecart(index,:) < dConcMax ...
               | t>tf
           OnContinue = false;
       end
       disp([' Temps : ' num2str(t/3600/24) 'j -' ...
             ' Numero de sauvegarde : ' num2str(index)  ...
             ' - Concentration Totale : ' num2str(sum(C*dx)) ...
             ' - Ecart : ' num2str(Ecart(index,:))])
             %' - Vitesse : ' num2str(u)])
       TempsConc(index)=t;

       figure(2)
       AffichageConcentration

       pause(0.01)
   end
end

Ccalc = C_analytical(mean(Ws), mean(Ks), z_, sum(C_init*dx), L);
disp(['MSE analytical // model : ', num2str( MSE(C,Ccalc) ), ' MP.m^-3'])
plot(Ccalc, -z_, 'g')