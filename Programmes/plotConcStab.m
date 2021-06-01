% load data as zHistory
% change parameters if needed

L = 67.3456;
N = 50;
dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

dt = 10;
tf = 60*60*24*10;
dtStab = 30*60;

timeLine = 0:dt:tf;
nCase = dtStab/dt;

[meanConc1, stdConc1] = getMeanConc(zHistory(end-2*nCase:end-nCase), length(z_), dz);
[meanConc2, stdConc2] = getMeanConc(zHistory(end-nCase:end), length(z_), dz);

figure(1), clf, hold on,

p1 = plot(meanConc1, -z_, 'b');
plot(meanConc1+2*stdConc1, -z_, '--b', meanConc1-2*stdConc1, -z_, '--b');
p2 = plot(meanConc2, -z_, 'r');
plot(meanConc2+2*stdConc2, -z_, '--r', meanConc2-2*stdConc2, -z_, '--r');

legend([p1, p2],{'tf-45min','tf-15min'}, 'Location', 'best');

xlabel('Concentration (mps.m⁻¹)')
ylabel('Depth (m)')

hold off