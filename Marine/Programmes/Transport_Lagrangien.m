% Ws qui minimise l'Err ...
% (Csurf peut etre extrapole par C(-1m) )
% 
% On connaît :
%     dx
%     C_
%     Z_
%     dh
global dt

Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];

% Initilisation :
D=350e-6; %m : Diametre
ModeleHydro='2012RHOMA_arome_003.nc';
indNom=3;

SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)

tf= 100*86400; dtmax=0.01; 
Tdes=60; %dxMax=0; % Tdes : intervalle de temps entre les tests d'équilibre 
% dConcMax : seuil de delta de concentration à partir duquel on de considère à l'équilibre 
dh=0.15; % profondeur sur laquelle le filet prélève

clear Concentration err
L = 50; % Profondeur (arbitraire)
N=2000;  dx= L/N;  x=0:dx:L; % x : boundaries of the meshes
x_=(x(1:end-1)+x(2:end))/2; % milieu de chaque maille
Xmin=0;Xmax=L;Cmin=-1;Cmax=2;
        %z0=H0(I0,J0)*Sigma;
        z0=L*Sigma;
        %z0_=H0(I0,J0)*(Sigma(1:end-1)+Sigma(2:end))/2;
        z0_=L*(Sigma(1:end-1)+Sigma(2:end))/2;    

% conditions initiales 
CMes=[0.62 0.34 0.06 0.02 0]; % concentrations mesurées
ZMes=[1 10 15 40 L]; % profondeur de chaque mesure
C = interp1(ZMes,CMes,x(1:end-1)+dx/2,'pchip'); % interpolation sur x
C=max(0*C,C); 

% Initialisation des positions de chaque particule
np = 100;
pd = makedist('Normal');
r = random(pd, 1, np);
part = min(L, abs(r)*L/2);
disp(part)
    
%% 1) calculer le profil de concentration associé à Ws
row = 1000;

% Determiner rho eau
DensiteFevrierRhoma
Nu=interp1(z0,KZ_Fev10,-x_,'pchip');
%Ks = 0.01;
%Nu = ones(1,2000)*Ks;

InitialisationVitesseTransport

rop=1011.4;
S=rop./row;     D_=((g*(abs(S-1))/nuw^2).^(1/3))*D;
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
u=Ws; u(rop<row)=-Ws(rop<row);

u0_=max(u);Nu0_=max(Nu);
if (u0_~=0 && Nu0_~=0) 
   dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
elseif (u0_==0 && Nu0_~=0) 
   dt=dx*dx/(2*Nu0_)*0.5;
elseif (u0_~=0 && Nu0_==0)
   dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
else
   dt=dtmax;
end

index = max(1, cast(part(1,:)/dx, 'uint32'));
part = [part; u(index)];
part = [part ; Nu(index)];

% part = [x1  x2 ... xn ; 
%         u1  u2 ... un ;
%        Nu1 Nu2 ... Nu3]

t=0; OnContinue=true;
while OnContinue
   t=t+dt;
    
   index = max(1, cast(part(1,:)/dx, 'uint32'));
   part(2,:) = u(index);
   part(3,:) = Nu(index);

   temp_part = part; % part(t-1)

   part(1,:) = arrayfun(@(i) Step_Lagrangien(part(2,i), part(3,i), part(1,i)), 1:size(part,2));
   

   if (mod(t,Tdes)<=dt/2 || Tdes-mod(t,Tdes)<=dt/2 )

       Ecart = abs(part(1,:) - temp_part(1,:));
       
       if (t>tf | part(1,:)<1E-9)
           OnContinue = false;
       end
       
       %disp(x_part);
       %disp(Ecart);
       %disp([' Temps : ' num2str(t/3600/24) 'j -' ...
       %      ' Compartiment : ' num2str(index)  ...
       %      ' - x : ' num2str(part(1,:)) ...
       %      ' - Ecart : ' num2str(Ecart)])
       disp(part)
   end
end

