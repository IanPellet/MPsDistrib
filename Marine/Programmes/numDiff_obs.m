%N_tested = [10 50 100 200 500 1000]; 
%N_tested = 1000:100:2000;
N_tested = 100:50:1000;
%[C_obt, z] = arrayfun(@(n) Transport_Eulerian(n), N_tested, 'UniformOutput', false);
leg = arrayfun(@(i) num2str(i), N_tested, 'UniformOutput', false);

L = 50;
N_max = N_tested(end);
dz_max= L/N_max;  
z_max=0:dz_max:L;

mse_list = ones(size(N_tested));

figure(1), clf;
hold on
i = 0;
%for i = 1:size(C_obt,2)
for N = N_tested
    %Ci = cell2mat(C_obt(i));
    %zi = cell2mat(z(i));
    %plot(Ci, -zi)
    i = i+1;
    
    [Cn, zn, MSEn] = Transport_Eulerian(N, L);
    Cni = interp1(zn, Cn, z_max, 'pchip');
    
    mse_list(i) = MSEn;
    
    plot(Cni, -z_max);
    
end
legend(leg)
xlabel('Concentration (MPs.m-3)')
ylabel('Profondeur (m)')
hold off

figure(2), clf;
semilogy(N_tested, mse_list)
xlabel('Number of meshes')
ylabel('Error')