dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  
nProfile = length(meanConc);
figure(1), clf,
hold on
for i = 100:10:nProfile
plot(meanConc{i},-z_,'DisplayName', [num2str((i*30-15)/60) ' h'])
end
legend('Location', 'best')
hold off