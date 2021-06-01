% Lagrangian transportation model 
% 
% 

global dt L
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
tf= 0.5*86400; % maximum simulation time
dt_max=0.01; % maximun time interval
dt_test = 0.1; % time interval between equilirium tests

%% Water column parameters
Lon0= 5.29;Lat0=43.24; %point Somlit
[I0,J0] = ReperePoint(Lon,Lat,Lon0,Lat0); % indices corresponding to the location
%L = H0(I0,J0); % depth
L = 50;
N=50;  dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh
z_min=0;z_max=L;C_min=-1;C_max=2;
z0=L*Sigma;
z0_=L*(Sigma(1:end-1)+Sigma(2:end))/2;    

%% Initial concentrations
dh=0.15; % depth of the net
CMes=[0.62 0.34 0.06 0.02 0]; % measured concentrations
ZMes=[1 10 15 40 L]; % Depth of each measure
C = interp1(ZMes,CMes,z(1:end-1)+dz/2,'pchip'); % interpolation on z
C=max(0*C,C); % negative values set to 0

%% Particules initialisation
D=350e-6; % Particle diametre
rop=1010.5; % Particle density
N_part0 = 50; % Number of particle wanted
P_part = C*dz/sum(C*dz); % percentage of particles in each mesh
n = round(N_part0*P_part); % nbr of part per mesh
N_part = sum(n); % actual nbr of part initialised
disp(['Total nbr of particles : ', num2str(N_part)])

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

Nmes = sum(C*dz);
figure(1),clf,
hold on
plot(CMes*dz*N_part/Nmes,-ZMes,'og', C*dz*N_part/Nmes,-z_,'r', n, -z_, 'b'),
title("Initial conditions", ["N_p_a_r_t = ",num2str(N_part)]),
legend(["Measured", "Interpolated", "Model init"])
hold off


%% Speed initialisation
DensiteFevrierRhoma 
%[K,dK] = wcp_interpolation(z0,KZ_Fev10,-z_); % Diffusivity
K_val = 0.001;
K = ones(size(z))*K_val;
dK = zeros(size(z));
InitialisationVitesseTransport
S=rop./row;     
D_=((g*(abs(S-1))/nuw^2).^(1/3))*D;
Ws=eval(['Vitesse' cell2mat(Nom(indNom)) '(D,S,D_);']);
u=Ws; u(rop<row)=-Ws(rop<row);


%% Particules matrix definition
% part = [z1 z2 ... zn ; 
%         u1 u2 ... un ;
%         K1 K2 ... Kn ;
%        dK1 dK2...dK3]
index = max(1, cast(z_part/dz, 'uint32'));
u_part = u(index);
K_part = K(index);
dK_part = dK(index);
part = [z_part ; u_part ; K_part ; dK_part];

%% Setting dt
u0_=max(u);Nu0_=max(K);
if (u0_~=0 && Nu0_~=0) 
   dt=min(dz/abs(u0_)*0.5,dz*dz/(2*Nu0_)*0.5); 
elseif (u0_==0 && Nu0_~=0) 
   dt=dz*dz/(2*Nu0_)*0.5;
elseif (u0_~=0 && Nu0_==0)
   dt=dz/abs(u0_)*0.5; 
else
   dt=dt_max;
end


%% Simulation
%figure(2), clf
h_init = histogram(part(1,:), "BinEdges", z, 'Visible', 'off').Values;
C_init = h_init/dz;
C_history = [C_init];
z_past = part(1,:);
t=0; OnContinue=true;

figure(4)
while OnContinue
    
    % Time update
    t=t+dt;
    
    % Particules update
    temp_part = part; % part(t-1)
    part(1,:) = Step_Lagrangien(part(1,:), part(2,:), part(3,:), part(4,:));
    %part(1,:) = step(part(1,:), part(2,:), part(3,:), part(4,:), alpha, dz);
    index = max(1, cast(part(1,:)/dz, 'uint32'));
    part(2,:) = u(index);
    part(3,:) = K(index);
    part(4,:) = dK(index);

    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 ) 
        % Equilibrium test
        z_present = part(1,:);
        % Computation of the concentration of MPs in each mesh
        %figure(2), clf
        h_past = histogram(z_past, "BinEdges", z, 'Visible', 'off').Values;
        h_present = histogram(z_present, "BinEdges", z, 'Visible', 'off').Values;
        C_past = h_past/dz;
        C_present = h_present/dz;
        C_history = [C_history ; C_present];
                
        dC = max(abs(C_present - C_past)/dt_test);
        

        if (t>tf)
           OnContinue = false;
        end

        disp([' Temps : ' num2str(t/3600/24) 'j -' ...
              ' - Ecart : ' num2str(dC)])
        

        %figure(3)
        h = histogram(part(1,:), "BinEdges", z, 'Visible', 'off' );
        C = h.Values/dz;

        %figure(2)
        %clf, hold on,
        %plot(C, -z_)
        %figure(4)
        plot(-part(1,:),'x')
        ylim([-L 0])
        grid on
        inv = z(end:-1:1);
        yticks(-inv)
        %plot(C/dz,-z_,'r',CMes/dz,-ZMes,'og', n, -z_, 'b')
        %plot(CMes,-ZMes,'og')
        pause(dt_test)
        
        z_past = z_present;
    end
end

%% Results
figure(5), clf, 

C = C_history(end,:);

Ccalc = C_analytical(mean(Ws), K_val, z_, N_part, L);
disp(sum(Ccalc*dz))
plot(Ccalc*dz, -z_, C*dz, -z_)

mse = MSE(C,Ccalc) ;
disp(['MSE analytical // model : ', num2str(sqrt(mse)), ' MPs.m⁻¹'])
