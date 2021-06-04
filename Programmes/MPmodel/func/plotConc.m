dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  
nProfile = length(meanConc);

nih = 2*24;
figure(1), clf,
hold on
for i = 1:nih:nProfile
plot(meanConc{i},-z_,'DisplayName', [num2str((i*30-15)/60) ' h'])
end
legend('Location', 'best')
hold off
xlabel('Concentration (mps.m⁻³)')
ylabel('Depth (m)')
title('Concentration profiles each day of a 10 days simulation')