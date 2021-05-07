% D_test = [125
% 225
% 375
% 750
% 1250
% 1750
% 2750
% 3750
% 4750
% ]*1e-6;
D_test = 350e-6;
RhoP = 1010.05;

resE = struct('C', [], 'z', [], 'D', 0);
resL = struct('C', [], 'z', [], 'D', 0);

type = 0;
nPart = 50e3;
tf = 60*60*24*2;
dt_test = 60*60*2;
L = 50;
N = 500;

for i=1:length(D_test)
% for i=1:length(RhoP_test)
    D = D_test(i);
    [resE(i).C, resE(i).z] = Transport_Eulerian(D, RhoP, N);
    [resL(i).z, resL(i).C, ~] = varMP_model(D, RhoP, type, nPart, tf, dt_test, 0, 0, 0, 0, L, N, '../');
    resE(i).D = D;
    resL(i).D = D;
end, clear i,

a = [resE(1).z];
b = [resL(1).z];
dz_E = a(2)-a(1);
dz_L = b(2)-b(1);

figure(2),clf,
hold on 
for i=1:length(D_test)
    alpha = (sum(resL(i).C)*dz_L)/(sum(resE(i).C)*dz_E);
%     plot([resE(i).C], -resE(i).z, 'DisplayName', ['Eul. D = ' num2str(resE(i).D*1e6) ' µm']);
%     plot([resL(i).C], -resL(i).z, 'DisplayName', ['Lag. D = ' num2str(resL(i).D*1e6) ' µm']);
    plot([resE(i).C]*alpha, -resE(i).z, 'DisplayName', ['Eulerian model']);
    plot([resL(i).C], -resL(i).z, 'DisplayName', ['Lagrangian model']);
end, clear i,
legend('Location', 'best')
xlabel('Concentration (mps.m⁻³)')
ylabel('Depth (m)')
% xlim([0 100])
ylim([-50 0])
title(['D = ' num2str(D_test(1)*1e6) 'µm -- rho = ' num2str(rhoP)...
    'kg.m⁻³ -- nPart = ' num2str(nPart) ' -- N = ' num2str(N)] )
hold off