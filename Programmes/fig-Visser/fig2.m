% Plot particles position concentration on the column during the simulation

Npart = 5000;
H = 50;
N = 50;
tf = 1e5;
dt_test = 60*60;

[TimeSimu,PartPos,z,z_,K] = RunSimu(Npart,H,N,tf,dt_test);

density = zeros(length(z_),length(TimeSimu));
for t = 1:length(TimeSimu)
    density(:,t) =  histogram(PartPos(t,:), "BinEdges", z, 'Visible', 'off').Values;
end

figure(1), clf,
ax1 = subplot(1,2,1);
plot(ax1,K,-z_)
xlabel('Diffusivity (m².s⁻¹)');
ylabel('Depth (m)');

ax2 = subplot(1,2,2);
pcolor(ax2, TimeSimu,-z_,density);
c = colorbar;
c.Label.String = 'Number of particles';
xlabel('Time (s)');
ylabel('Depth (m)');
