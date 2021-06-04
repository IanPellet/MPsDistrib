clear

load('/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/Data/2021-6-2-20-39-51_234344-zHist864000.mat')
load('/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/2021-6-2-20-39-51_234344.mat')

dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

nProfile = length(zHistory);
tinit = tf - nProfile*dt;
nCase = dtAvgC/dt;

[meanConc, stdConc, cHistory] = getMeanConc(zHistory(nProfile-nCase+1:nProfile), length(z_), dz);

figure(1), clf,
hold on
for i = (tf-30*60-tinit)/dt:nProfile
plot(cHistory(i,:),-z_,'DisplayName', [num2str(tinit+i*dt) ' s'])
end
% legend('Location', 'best')
hold off
xlabel('Concentration (mps.m⁻³)')
ylabel('Depth (m)')
title('Concentration profiles over the last 30min of a 10 days simulation', ['dt = ' num2str(dt) ' s']) 

figure(2), clf,
hold on
plot(meanConc, -z_, 'b')
plot(meanConc+2*stdConc, -z_, '--b')
plot(meanConc-2*stdConc, -z_, '--b')
hold off
xlabel('Concentration (mps.m⁻³)')
ylabel('Depth (m)')
title('Time average of concentration profiles over 30min', 'C_a_v_g +/- 2std') 

