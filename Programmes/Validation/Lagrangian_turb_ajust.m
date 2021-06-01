% Values taken by alpha
min_a = 0.5;
max_a = 1.5;
n_a = 50;
alpha = linspace(min_a, max_a, n_a);

% Values taken by Kz
ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro)
L = 50; % depth
N=100;  dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh
DensiteFevrierRhoma 
[K,~] = wcp_interpolation(z0,KZ_Fev10,-z_); % Diffusivity

min_K = min(K);
%max_K = max(K);
max_K = 0.04;
n_K = 10; % number of k tested 
K_test = linspace(min_K, max_K, n_K);

% Initialisation of the error values list
mse_list = ones(size(K_test,2),size(alpha,2));
min_alpha = zeros(1,size(K_test,2));
min_error = zeros(1,size(K_test,2));

% file name format
formatOut = 'mmdd';
f_0 = ['../../Ian/Results/MSE_alpha/',datestr(now,formatOut),'_mse-a_a',num2str(min_a),'-',...
        num2str(max_a),'_Np'];

i_f = 0;
f_test = fullfile([f_0,'*',num2str(n_K),'_0','.*']);
while exist(f_test, 'file')
    i_f = i_f+1;
    f_test = fullfile([f_0,'*',num2str(n_K),'_',num2str(i_f),'.*']);
end

f_end = ['-',num2str(n_K),'_',num2str(i_f)];


i_K = 0;
for k = K_test
    i_K = i_K+1;
    fprintf(['\n\n--------------------- Ks = ' num2str(k) ' ---------------------\n'])
    
    i_a = 0;
    for a = alpha
        i_a = i_a+1;
        fprintf(['\n\n--------------------- Alpha = ' num2str(a) ' ---------------------\n\n'])
        figure(1); clf;
        [Ca, MSEa, z, Np] = transport(a, k);
        mse_list(i_K,i_a) = MSEa;
    end
   
    f = figure(i_K+1); clf;
    semilogy(alpha,sqrt(mse_list(i_K,:)))
    title("Error between Lagrangian model and Analytical solution",...
        ['K_s = ',num2str(k),' m².s⁻¹, N_p_a_r_t = ',num2str(Np)])
    xlabel('Alpha')
    ylabel('Error (MPs.m⁻¹)')
        
    % Export figure as eps file
    f_name = [f_0,num2str(Np),'_K',num2str(i_K),f_end,'.eps'];
    removeToolbarExplorationButtons(f);
    exportgraphics(f,f_name,'ContentType','vector');
    
    min_a_k = min(alpha(mse_list(i_K,:) == min(mse_list(i_K,:))));
    disp(min_a_k)
    min_alpha(i_K) = min_a_k;
    min_error(i_K) = min(mse_list(i_K,:));
end

min_a_Ks = [K_test ; min_alpha ; min_error];
fvar_name = [f_0,num2str(Np),'_K',num2str(i_K),f_end,'.mat'];
save(fvar_name, 'min_a_Ks')
disp(min_a_Ks)

%% STEP 
function [new_z] = step(z, u, Kz, dKz, alpha, dz)
%STEP One step foward in time for the Lagrangian transport model
%   The new value of z after dt is computed and returned by this function.
    global dt L

    d = 1;
    pd = makedist('Normal', 'mu', 0, 'sigma', d);
    R = random(pd,size(z));
    new_z = z + dt*u + dKz*dt/dz + alpha*R.*(2/d*Kz*dt).^(1/2);
    %new_x = max(0,new_x);
    %new_x = min(L, new_x);

    % Reflective boundary conditions
    new_z(new_z < 0) = -new_z(new_z < 0);
    new_z(new_z > L) = L-new_z(new_z > L) + L;
end


%% MODEL
function [C, mse, z_, N_part] = transport(alpha, K_val)
%TRANSPORT Lagrangian transportation model 

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
    tf= 0.2*86400; % maximum simulation time
    dt_max=0.01; % maximun time interval
    dt_test = 60*60; % time interval between equilirium tests

    %% Water column parameters
    %Lon0= 5.29;Lat0=43.24; %point Somlit
    %[I0,J0] = ReperePoint(Lon,Lat,Lon0,Lat0); % indices corresponding to the location
    L = 50; % depth
    N=L;  
    dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
    z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh
    z0=L*Sigma; 

    %% Initial concentrations
    CMes=[0.62 0.34 0.06 0.02 0]; % measured concentrations
    ZMes=[1 10 15 40 L]; % Depth of each measure
    C = interp1(ZMes,CMes,z(1:end-1)+dz/2,'pchip'); % interpolation on z
    C=max(0*C,C); % negative values set to 0

    %% Particules initialisation
    D=350e-6; %m : Diametre
    rop=1010.5;
    N_part0 = 1000;
    P_part = C*dz/sum(C*dz);
    n = round(N_part0*P_part);
    %n = round(C*N_part0/sum(C));
    N_part = sum(n);
    %n = round(C/dz); % number of particules per mesh
    %N_part = sum(n); % Total number of part in the water column
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


    %% Speed initialisation
    DensiteFevrierRhoma 
    %[K,dK] = wcp_interpolation(z0,KZ_Fev10,-z_); % Diffusivity
    K = ones(size(z_))*K_val;
    dK = zeros(size(z_));
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
    h_init = histogram(part(1,:), "BinEdges", z).Values;
    C_init = h_init/dz;
    C_history = [C_init];
    z_past = part(1,:);
    t=0; OnContinue=true;
    while OnContinue

        % Time update
        t=t+dt;

        % Particules update
        temp_part = part; % part(t-1)
        part(1,:) = step(part(1,:), part(2,:), part(3,:), part(4,:), alpha, dz);
        index = max(1, cast(part(1,:)/dz, 'uint32'));
        part(2,:) = u(index);
        part(3,:) = K(index);
        part(4,:) = dK(index);

        if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 )

            % Save state to history



            % Equilibrium test
            z_present = part(1,:);
            % Computation of the concentration of MPs in each mesh
            h_past = histogram(z_past, "BinEdges", z).Values;
            h_present = histogram(z_present, "BinEdges", z).Values;
            C_past = h_past/dz;
            C_present = h_present/dz;
            C_history = [C_history ; C_present];

            dC = max(abs(C_present - C_past)/dt_test);


            if (t>tf || dC == 0)
               OnContinue = false;
            end

            %disp([' Temps : ' num2str(t/3600/24) 'j -' ...
            %      ' - Ecart : ' num2str(dC)])

            z_past = z_present;
        end
    end

    C = C_history(end,:);

    Ccalc = C_analytical(mean(Ws), mean(K_val), z_, N_part, L);
    disp(sum(Ccalc*dz))
    plot(Ccalc, -z_, C, -z_)
    title(["alpha = ",alpha])
    mse = MSE(C,Ccalc) ;
    disp(['MSE analytical // model : ', num2str(sqrt(mse)), ' MPs.m⁻¹'])
end