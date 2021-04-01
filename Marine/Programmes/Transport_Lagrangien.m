% Lagrangian transportation model 
% 
% 

global dt

%% Speed computation formulas
Nom=[... 
    ;{'Nielsen'}...%Nielsen (1992)
    ;{'Soulsby'}...%Soulsby (1997)
    ;{'Ahrens'}...%Ahrens (2000)
   ];

D=350e-6; %m : Diametre
ModeleHydro='2012RHOMA_arome_003.nc';
indNom=3;

SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)

%% Equilibrium test parameters
tf= 100*86400; % maximum simulation time
dt_max=0.01; % maximun time interval
dt_test = 60*30; % time interval between equilirium tests
dC_min = 5E-5; % C(t+dt)-C(t) threshold for the system to be considered at equilibrium

clear Concentration err

%% Water column parameters
dh=0.15; % depth of the net
L = 50; % depth
N=2000;  dx= L/N;  x=0:dx:L; % x : boundaries of the meshes
x_=(x(1:end-1)+x(2:end))/2; % middle of each mesh
Xmin=0;Xmax=L;Cmin=-1;Cmax=2;
        %z0=H0(I0,J0)*Sigma;
        z0=L*Sigma;
        %z0_=H0(I0,J0)*(Sigma(1:end-1)+Sigma(2:end))/2;
        z0_=L*(Sigma(1:end-1)+Sigma(2:end))/2;    

%% Initial concentrations
CMes=[0.62 0.34 0.06 0.02 0]; % measured concentrations
ZMes=[1 10 15 40 L]; % Depth of each measure
C = interp1(ZMes,CMes,x(1:end-1)+dx/2,'pchip'); % interpolation on x
C=max(0*C,C); % negative values set to 0

%% Particules position initialisation
n = round(C/dx/2); % number of particules per mesh
N_part = sum(n); % Total number of part in the water column
x_part = ones(1, N_part); % Position of each part ; space allocation
i_part = 0; % Part index
for i = 1:N % Mesh index
    % Each particle is initalised at a random depth of its mesh
    % Probability distribution : Uniform distribution in the mesh
    pd = makedist('Uniform','lower',x(i),'upper',x(i+1)); 
    r = random(pd, 1, n(i)); 
    temp_j = i_part+1;
    i_part = i_part+n(i);
    x_part(temp_j:i_part) = r; % x_part is filled with the position of each part
end

%h = histogram(part(1,:), "BinEdges", x);
%C = h.Values*dx;

figure(1),clf,plot(C/dx,-x_,'r',CMes/dx,-ZMes,'og', n, -x_, 'b');



% Determiner rho eau
row = 1000;
DensiteFevrierRhoma
Nu=interp1(z0,KZ_Fev10,-x_,'pchip');

InitialisationVitesseTransport

rop=1011.4;
S=rop./row;     D_=((g*(abs(S-1))/nuw^2).^(1/3))*D;
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
u=Ws; u(rop<row)=-Ws(rop<row);


%% Particules matrix definition
% part = [x1  x2 ... xn ; 
%         u1  u2 ... un ;
%        Nu1 Nu2 ... Nu3]
index = max(1, cast(x_part/dx, 'uint32'));
u_part = u(index);
Nu_part = Nu(index);
part = [x_part ; u_part ; Nu_part];

%% Setting dt
u0_=max(u);Nu0_=max(Nu);
if (u0_~=0 && Nu0_~=0) 
   dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
elseif (u0_==0 && Nu0_~=0) 
   dt=dx*dx/(2*Nu0_)*0.5;
elseif (u0_~=0 && Nu0_==0)
   dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
else
   dt=dt_max;
end

%% Simulation
C_history = [];
x_past = part(1,:);
t=0; OnContinue=true;
while OnContinue
    
    % Time update
    t=t+dt;
    
    % Particules update
    temp_part = part; % part(t-1)
    part(1,:) = Step_Lagrangien(part(1,:), part(2,:), part(3,:));
    index = max(1, cast(part(1,:)/dx, 'uint32'));
    part(2,:) = u(index);
    part(3,:) = Nu(index);

    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 )
        
        % Save state to history
        C_history = [C_history ; part(1,:)];
        
        
        % Equilibrium test
        x_present = part(1,:);
        % Computation of the concentration of MPs in each mesh
        h_past = histogram(x_past, "BinEdges", x).Values;
        h_present = histogram(x_present, "BinEdges", x).Values;
        C_past = h_past*dx;
        C_present = h_present*dx;
                
        dC = max(abs(C_present - C_past)/dt_test);
        

        if (t>tf || dC < dC_min)
           OnContinue = false;
        end

        disp([' Temps : ' num2str(t/3600/24) 'j -' ...
              ' - Ecart : ' num2str(dC)])
        

        figure(3)
        h = histogram(part(1,:), "BinEdges", x);
        C = h.Values*dx;

        figure(2)
        clf, hold on,
        plot(C/dx, -x_)
        figure(4)
        clf
        plot(-part(1,:))
        %plot(C/dx,-x_,'r',CMes/dx,-ZMes,'og', n, -x_, 'b')
        %plot(CMes,-ZMes,'og')
        pause(0.01)
        
        x_past = x_present;
    end
end

%% Results

