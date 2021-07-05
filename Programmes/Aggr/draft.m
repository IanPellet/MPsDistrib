H = 50; % Depth of the column (m)
dz = 1e-3;

% x = 0:1e-2:1; % (m)
% y = 0:1e-2:1; % (m)
z = 0:dz:H; % (m)

Npart = 1e3;
X = rand(Npart,1);
Y = rand(Npart,1);
Z = rand(Npart,1).*H;
pd = makedist('Normal', 'mu', 350e-6, 'sigma', 50e-6);
sizeP = random(pd,Npart,1);
clear pd,

rhop = 1025;

mp = getMPlist(Npart, sizeP, rhop, Rhow(z, H), 0);

tf = 60*60*24*1;
dt_test = 60*60*2;
tsave = tf;

[K, dK] = turbulence(z, dz);

[Final,dt] = Aggr_simulator(mp, X, Y, Z, K, dK, H, dz, tf, dt_test, 0);


% t = 0:dt:tf;
% % scatter3([Final{:,1}], [Final{:,2}], -[Final{:,3}], 5, t)
% xF = cell2mat(Final(:,1));
% yF = cell2mat(Final(:,2));
% zF = cell2mat(Final(:,3));
% 
% for i = 1:Npart
%     plot3(xF(:,i), yF(:,i), -zF(:,i))
%     hold on 
% end
% hold off
% zlim([-H 0])
% ylim([0 1])
% xlim([0 1])


function [rho] = Rhow(z, H)
    Rm = 1.028505199810085e+03;
    RM = 1.028660884343125e+03;
    
    rho = Rm + z.*(RM-Rm)/H;
end