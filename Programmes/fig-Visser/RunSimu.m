function [TimeSimu,PartPos,z,z_,K] = RunSimu(eq,Npart,H,N,tf,dt_test)
%%RUNSIMU Run simulation of the Lagrandian diffusivity model

if nargin < 6
    dt_test = 60*60*2; % test time interval (s)
end
if nargin < 5 
    tf = 1e5; % simulation time (s)
end
if nargin < 4
    N = 500;
end
if nargin < 3
    H = 50;
end
if nargin < 2
    Npart = 5000;
end
if nargin < 1
    eq = 'eq3_NaiveRandomWalk';
end

%% Water column
dz= H/N;  z=0:dz:H; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh

[K, dK] = Diffusivity(z,z_,dz);
% K = 0.5./z_;
% dK = ones(size(K))*0.5/dz;

%% Particles
z_part = linspace(0, H, Npart);
K_part = zeros(size(z_part));
dK_part = zeros(size(z_part));
Part = [z_part ; K_part ; dK_part];
Part = UpdatePart(Part,K,dK,dz);

%% Time parameters
dt = 10;
ddK = diff(dK)/dz;
dt = min(dt, abs(min(1./ddK)/10)); % condition dt<<min(1/ddK) 
% dt = min(dt, dz/max(abs(u))); % condition dt < dz/max|u|

if rem(tf,dt_test) == 0
    nStep = tf/dt_test;
else
    nStep = fix(tf/dt_test)+1;
end

t = 0;
step = 1;

%% Output parameters
PartPos = -ones(nStep,Npart);
PartPos(step,:) = Part(1,:);

TimeSimu = ones(nStep,1);
TimeSimu(step) = t;

%% Simulation
OnContinue=true;
while OnContinue
    % Time update
    t=t+dt;
        
    % Particules update
    zi = Part(1,:);
    Kzi = Part(2,:);
    dKzi = Part(3,:);
%     newz = eq3_NaiveRandomWalk(zi,Kzi,dKzi,K,dt,dz);
    newz = eval([eq '(zi,Kzi,dKzi,K,dt,dz)']);
    Part(1,:) = newz;
    % Boundary conditions
    Part(1,newz < 0) = -newz(newz < 0);
    Part(1,newz > H) = H - newz(newz > H) + H;
    
    Part = UpdatePart(Part,K,dK,dz);
    
    
    if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 ) % test
        step = step+1;
        TimeSimu(step) = t;
        PartPos(step,:) = Part(1,:);
        
        if (t>=tf)
           OnContinue = false;
        end
        
        dC = DeltaConcentration(PartPos(step-1,:), PartPos(step,:), z, dz,dt);
        disp([' Temps : ' num2str(t/tf*100) '% - Ecart : ' ...
            num2str(dC*100), '%'])

    end
end
end