% Plot particles position concentration on the column during the simulation

ModeleHydro='2012RHOMA_arome_003.nc';
SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
load(SauvegardeModeleHydro,'Lon','Lat','H0')

% Npart = 5000;
Lon0= 5.29;Lat0=43.24; %point Somlit
[I0,J0] = ReperePoint(Lon,Lat,Lon0,Lat0); % indices corresponding to the location
H = H0(I0,J0); % depth
% H = 55;
N = fix(H);
Npart = N*100;
tf = 1e5;
dt_test = 60*60;
eq = 'eq6_CorrectedRandomWalk';

[TimeSimu,PartPos,z,z_,K,dK] = RunSimu(eq,Npart,H,N,tf,dt_test);

density = zeros(length(z_),length(TimeSimu));
for t = 1:length(TimeSimu)
    density(:,t) =  histogram(PartPos(t,:), "BinEdges", z, 'Visible', 'off').Values;
end

Ref = density(:,1);
err = zeros(size(TimeSimu));
for t = 1:length(TimeSimu)
    err(t) =  mean(abs((density(:,t)-Ref)./Ref))*100;
end

nFig = 4;
eq2 = split(eq,'_');
eq3 = [cell2mat(eq2(1)) ' ' cell2mat(eq2(2))];
% ttl = [eq3 ' with ' num2str(Npart) ' particles'];
% title(ttl)
figure('WindowState','maximized', 'name', [cell2mat(eq2(1)) '-' num2str(Npart) '-ptSomlit']), clf,

ax1 = subplot(1,nFig,1);
% plot(ax1,K,-z_)
% xlabel('Diffusivity (m².s⁻¹)');
% ylabel('Depth (m)');
plot(ax1,TimeSimu,err)
ylabel('Mean Error (%)');
xlabel('Time (s)');

ax2 = subplot(1,nFig,2);
pcolor(ax2, TimeSimu,-z_,density);
c = colorbar;
c.Label.String = 'Number of particles';
xlabel('Time (s)');
ylabel('Depth (m)');

ax3 = subplot(1,nFig,3);
plot(ax3,(density(:,end)-Ref)./Ref*100,-z_)
ylim([-H,0]);
xlabel('Error (%)');
ylabel('Depth (m)');

ax4 = subplot(1,nFig,4);
plot(ax4,K,-z_)
ylim([-H,0]);
xlabel('Diffusivity (m².s⁻¹)');
ylabel('Depth (m)');

