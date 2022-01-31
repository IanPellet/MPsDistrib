clear,
figPath = '/media/ian/Transcend/MPsDistrib/Results/';

BeaufortForce = 4;
% v=3.B^{3/2} si v  est en km/h, Wikipedia
wind_speed = 3*BeaufortForce.^(3/2);
H = 70;
tf = 60*60*24*3; % Trial time
P = 15000; % Number of particles
dt_test = 60*60;
partZinit = linspace(0, H, P);

Psizes = ones(1,P)*400e-6;

% allocate memory to store particles
part_list(P) = Particle; % array of Particle objects
rhoMP = 1000;

load(['/media/ian/Transcend/MPsDistrib/Data/Turbulence/DataTurbForce' num2str(BeaufortForce) '.mat'])

% Sea water density profile
date = datetime(2020,03,18);
[~,Row_day,~,z__day] = KsSalTemp2020(wind_speed, date);
rhof = interp1(-z__day,Row_day,z,'pchip'); % density of sea water [kg.m⁻³]
clear KZ_day Row_day z_day z__day,


%% create list of particles
% Fill the array
for i = 1:(P)
    part_list(i) = Particle(Psizes(i), 1025, rhoMP, 0, i, min(Psizes(i))); % size and densities only used if ws computed
end, clear i,

[zFinal, ~, ~] = Part_Simulator(part_list, partZinit, K, dK, H, dz, rhof, tf, dt_test, 60, 1e-8, true, false, false, true, 60*60*2);
[ConcMean, ConcSTD, hConc] = getMeanConc(zFinal, length(z_), dz);

    
    %% Plot
f1 = figure(1); clf,
for iS = 1:length(zFinal)
    plot(hConc(iS,:),-z_)
    hold on    
end, clear iS,
xlabel('Concentration (particles.m⁻¹)')
ylabel('Depth (m)')

f1Name = 'ConcProf2h';
savefig(f1, [figPath f1Name '.fig']);
exportgraphics(f1, [figPath f1Name '.png']);

f2 = figure(2); clf,
plot(ConcMean,-z_,'k','DisplayName', 'Time average', 'lineWidth', 1.5)
hold on
plot(ConcMean+ConcSTD,-z_,'--','DisplayName', 'Average + std','color','k');
plot(ConcMean-ConcSTD,-z_,'--','DisplayName', 'Average - std','color','k');        
xlabel('Concentration (particles.m⁻¹)')
ylabel('Depth (m)')
legend('Location','best')
xlim([0 max(ConcMean+ConcSTD)])
    
f2Name = 'ConcProf2hAvgSTD';
savefig(f2, [figPath f2Name '.fig']);
exportgraphics(f2, [figPath f2Name '.png']);