% Lagrangian transportation model 
% 
% 

global dt
clear Concentration err

%% Speed computation formulas
Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];
indNom=3;

%% Load hydrodinamic model
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)

%% Equilibrium test parameters
tf= 100*86400; % maximum simulation time
dt_max=0.01; % maximun time interval
dt_test = 60*30; % time interval between equilirium tests
dC_min = 5E-5; % C(t+dt)-C(t) threshold for the system to be considered at equilibrium

%% Water column parameters
Lon0= 5.29;Lat0=43.24; %point Somlit
[I0,J0] = ReperePoint(Lon,Lat,Lon0,Lat0); % indices corresponding to the location
L = H0(I0,J0); % depth
dh=0.15; % depth of the net
N=2000;  dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh
z_min=0;z_max=L;C_min=-1;C_max=2;
z0=L*Sigma;
z0_=L*(Sigma(1:end-1)+Sigma(2:end))/2;    

%% Initial concentrations
CMes=[0.62 0.34 0.06 0.02 0]; % measured concentrations
ZMes=[1 10 15 40 L]; % Depth of each measure
C = interp1(ZMes,CMes,z(1:end-1)+dz/2,'pchip'); % interpolation on z
C=max(0*C,C); % negative values set to 0

%% Particules initialisation
D=350e-6; %m : Diametre
rop=1011.4;
n = round(C/dz/2); % number of particules per mesh
N_part = sum(n); % Total number of part in the water column
z_part = ones(1, N_part); % Position of each part ; space allocation
i_part = 0; % Part index
for i = 1:N % Mesh index
    % Each particle is initalised at a random depth of its mesh
    % Probability distribution : Uniform distribution in the mesh
    pd = makedist('Uniform','lower',z(i),'upper',z(i+1)); 
    r = random(pd, 1, n(i)); 
    temp_j = i_part+1;
    i_part = i_part+n(i);
    z_part(temp_j:i_part) = r; % z_part is filled with the position of each part
end

%h = histogram(part(1,:), "BinEdges", z);
%C = h.Values*dz;

figure(1),clf,plot(C/dz,-z_,'r',CMes/dz,-ZMes,'og', n, -z_, 'b');



%% Speed initialisation
DensiteFevrierRhoma 
Nu=interp1(z0,KZ_Fev10,-z_,'pchip'); % Diffusivity
InitialisationVitesseTransport
S=rop./row;     D_=((g*(abs(S-1))/nuw^2).^(1/3))*D;
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
u=Ws; u(rop<row)=-Ws(rop<row);


%% Particules matrix definition
% part = [z1 z2 ... zn ; 
%         u1 u2 ... un ;
%         K1 K2 ... kn]
index = max(1, cast(z_part/dz, 'uint32'));
u_part = u(index);
Nu_part = Nu(index);
part = [z_part ; u_part ; Nu_part];

%% Setting dt
u0_=max(u);Nu0_=max(Nu);
if (u0_~=0 && Nu0_~=0) 
   dt=min(dz/abs(u0_)*0.5,dz*dz/(2*Nu0_)*0.5); 
elseif (u0_==0 && Nu0_~=0) 
   dt=dz*dz/(2*Nu0_)*0.5;
elseif (u0_~=0 && Nu0_==0)
   dt=min(dz/abs(u0_)*0.5,dz*dz/(2*Nu0_)*0.5); 
else
   dt=dt_max;
end

%% Simulation
C_history = [];
z_past = part(1,:);
t=0; OnContinue=true;
while OnContinue
    
    % Time update
    t=t+dt;
    
    % Particules update
    temp_part = part; % part(t-1)
    part(1,:) = Step_Lagrangien(part(1,:), part(2,:), part(3,:));
    index = max(1, cast(part(1,:)/dz, 'uint32'));
    part(2,:) = u(index);
    part(3,:) = Nu(index);

    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 )
        
        % Save state to history
        C_history = [C_history ; part(1,:)];
        
        
        % Equilibrium test
        z_present = part(1,:);
        % Computation of the concentration of MPs in each mesh
        h_past = histogram(z_past, "BinEdges", z).Values;
        h_present = histogram(z_present, "BinEdges", z).Values;
        C_past = h_past*dz;
        C_present = h_present*dz;
                
        dC = max(abs(C_present - C_past)/dt_test);
        

        if (t>tf || dC < dC_min)
           OnContinue = false;
        end

        disp([' Temps : ' num2str(t/3600/24) 'j -' ...
              ' - Ecart : ' num2str(dC)])
        

        figure(3)
        h = histogram(part(1,:), "BinEdges", z);
        C = h.Values*dz;

        figure(2)
        clf, hold on,
        plot(C/dz, -z_)
        figure(4)
        clf
        plot(-part(1,:))
        %plot(C/dz,-z_,'r',CMes/dz,-ZMes,'og', n, -z_, 'b')
        %plot(CMes,-ZMes,'og')
        pause(0.01)
        
        z_past = z_present;
    end
end

%% Results

