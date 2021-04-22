global nuw g row 
%global aW 

%%---Index des Vitesses et Legendes---%%

Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];

FacteurForme=[... 
    ;0....%Nielsen (1992)
    ;0....%Soulsby (1997)
    ;0....%Ahrens (2000)
   ];
FacteurForme=logical(FacteurForme);

NomLegende=[... 
           ;{'Nielsen'}...%Nielsen (1992)
           ;{'Soulsby (1997)'}...%Soulsby (1997)
           ;{'Ahrens (2000)'}...%Ahrens (2000)
          ];
      
Coul=[{'-b'},{'-c'},{'-m'},{'-r'},{'-k'},{'-g'},{'-y'}...
      ,{'-.b'},{'-.c'},{'-.m'},{'-.r'},{'-.k'},{'-.g'},{'-.y'}];

  
%%---Dimentionless grain size---%%
% D_=1:10;                D_3=D_.^3;


%%---Physical constants---%%

g=9.81 ; %m.s-1 (gravitational acceleration)
nuw=1.1*10^-6; %m2.s-1 (kinematic viscosity of sea water)
nuw2=nuw^2;

%%---EDM---%%
% Determiner rho eau
% DensiteFevrierRhoma
% row=1030;
row=interp1(-z__Fev10,Row_Fev10,z,'pchip'); %densite eau RHOMA10Fev
% row=interp1(-z__Fev03,Row_Fev03,x,'pchip'); %densite eau RHOMA03Fev
% row=interp1(-z__Avr16,Row_Avr16,x,'pchip'); %densite eau RHOMA16Avr

% %%---Particule---%%
% rop=1011.78; %sinking mps
% % rop=1011.69; %buoyant mps
% S=rop./row;
% 
% % d=D_./((g*(abs(S-1))/nuw^2).^(1/3)) ; 
% d=350e-6; %m
% D_=((g*(abs(S-1))/nuw^2).^(1/3))*d;
% D_3=D_.^3;



