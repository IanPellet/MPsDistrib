clear,
S = 1e-5;
H = 70;
N = 2*H;
%% INIT PARAMS
BeaufortForce = [0:4 8];
% BeaufortForce = 2;
% v=3.B^{3/2} si v  est en km/h, Wikipedia
wind_speed = 3*BeaufortForce.^(3/2);

tfInter = 60*60*24*3;
tf = 60*60*24*2; % simulation time
% addAggT = tf - 60*60*10;
saveLastSec = tf; 
dt_test = 60;
dt = 60;

%% create list of particles
rhoMP = 1000;
rhoB = 1100; %  doi:10.3354/meps175087

nMP = 200;
sizeMP = ones(nMP,1)*200e-6;
% zInitMP = ones(nMP,1)*20;

nAgg = 200;
sizeAgg = ones(nAgg,1)*1e-3; % doi.org/10.1038/srep00716
% zInitAgg = ones(nAgg,1)*30;

%% create list of MP
% allocate memory to store particles
mp_list(nMP) = Particle; % array of Particle objects
% Fill the array
for i = 1:(nMP)
    mp_list(i) = Particle(sizeMP(i), rhoB, rhoMP, 0, i, min(sizeMP));
end, clear i,
zMP = random('uniform',0,H,1,nMP);

% %% create list of Agrr
% % allocate memory to store particles
agg_list(nAgg) = Particle; % array of Particle objects
% Fill the array
for i = 1:(nAgg)
        agg_list(i) = Particle(sizeAgg(i), rhoB, rhoMP, 1, i, min(sizeMP));
end, clear i,
% zAgg = random('uniform',0,H,1,nAgg);

P = nMP+nAgg;
% partZinit = [zMP zAgg];

%% 
% part_list = [mp_list agg_list];
% clear mp_list agg_list sizeMP sizeAgg ,

for iB = [3 2 4:length(BeaufortForce)]
runID = strrep(join(split(num2str(clock)),'-'),'.','_');
load(['/media/ian/Transcend/MPsDistrib/Data/Turbulence/DataTurbForce' num2str(BeaufortForce(iB)) '.mat'])
H = 70;

date = datetime(2020,03,18);
[~,Row_day,~,z__day] = KsSalTemp2020(wind_speed(iB), date);
rhof = interp1(-z__day,Row_day,z,'pchip'); % density of sea water [kg.m⁻³]
clear Row_day z__day date u_result v_result K2,

[zInterMP, ~, MPinter,  ~, ~] = Part_Simulator2(mp_list, zMP, K, dK, H,...
    dz, rhof, tfInter, dt_test, dt, S, true, false, false, true, saveLastSec, 0,runID,false);
% [zInterAgg, ~, Agginter,  ~, ~] = Part_Simulator2(agg_list, zAgg, K, dK, H,...
%     dz, rhof, tfInter, dt_test, dt, S, true, false, false, true, saveLastSec, 0,runID,false);

pli = [MPinter agg_list];
part_listInter(P) = Particle;
for i = 1:(P)
    part_listInter(i) = Particle(pli(i).L, pli(i).Rhob, pli(i).RhoMP, pli(i).Fb, i, pli(i).L0);
end, clear i,
while isempty(zInterMP{end})
    zInterMP = zInterMP(1:end-1);
end
% while isempty(zInterAgg{end})
%     zInterAgg = zInterAgg(1:end-1);
% end
saveTo = '/media/ian/Transcend/MPsDistrib/Results/AggPlot/';
save([saveTo runID{:} '_simData.mat'])
% partZinter = [zInterMP{end} zInterAgg{end}];
partZinter = [zInterMP{end} random('uniform',0,40,1,nAgg)];
[zFinal, dtf, partFin, tauxAgg, ~, AggContent] = Part_Simulator2(part_listInter, partZinter, K, dK, H,...
    dz, rhof, tf, dt_test, dt, S, true, true, false, true, saveLastSec, 0,runID,true);
while isempty(zFinal{end})
    zFinal = zFinal(1:end-1);
    tauxAgg = tauxAgg(1:end-1);
end

save([saveTo runID{:} '_simData.mat'])

% %% Compute conc
% MPConcInter = GetConc(zInterMP, 1, nMP, N, dz);
% AggConcInter = GetConc(zInterAgg, 1, nAgg, N, dz);
% hConc = GetConc(zFinal, 1, nMP+nAgg, N, dz);
% MPConc = GetConc(zFinal, 1, nMP, N, dz);
% AggConc = GetConc(zFinal, nMP+1, nMP+nAgg, N, dz);
% 
% MeanContent = GetMeanContent(zFinal, nMP+1, nMP+nAgg, N, dz, AggContent);
% 
% save([saveTo runID{:} '_simData.mat'])
% 
% %% Plot
% % [allConcMat, allTauxMat, MPConcMat, AggConcMat] = plotAggPcolor(runID{:}, N, dz, nMP, nAgg);
% tminter = 0:dt:tfInter;
% timeline = linspace(0,tf,length(tauxAgg));
% 
% f1 = figure(1); clf,
% yyaxis left
% h = pcolor(timeline/60/60,-z_,hConc');
% % h = pcolor(tt'/60/60,-z,meancc');
% cb = colorbar;
% % caxis([0 max(max(meanC))])
% set(h, 'EdgeColor', 'none');
% xlabel('Time (h)')
% ylabel('Depth (m)')
% cb.Label.String = 'Concentration (particles.m⁻¹)';
% hold on
% yyaxis right
% plot(timeline/60/60,tauxAgg,'-')
% hold off
% title('All particles - Aggregation on')
% 
% 
% f2 = figure(2); clf,
% yyaxis left
% h = pcolor(timeline/60/60,-z_,MPConc');
% % h = pcolor(tt'/60/60,-z,meancc');
% cb = colorbar;
% % caxis([0 max(max(meanC))])
% set(h, 'EdgeColor', 'none');
% xlabel('Time (h)')
% ylabel('Depth (m)')
% cb.Label.String = 'Concentration (particles.m⁻¹)';
% hold on
% yyaxis right
% plot(timeline/60/60,tauxAgg,'-')
% hold off
% title('MP - Aggregation on')
% 
% 
% f3 = figure(3); clf,
% yyaxis left
% h = pcolor(timeline/60/60,-z_,AggConc');
% % h = pcolor(tt'/60/60,-z,meancc');
% cb = colorbar;
% % caxis([0 max(max(meanC))])
% set(h, 'EdgeColor', 'none');
% xlabel('Time (h)')
% ylabel('Depth (m)')
% cb.Label.String = 'Concentration (particles.m⁻¹)';
% hold on
% yyaxis right
% plot(timeline/60/60,tauxAgg,'-')
% hold off
% title('Agg - Aggregation on')
% 
% f4 = figure(4); clf,
% h = pcolor(flipud(MPConcInter'));
% cb = colorbar;
% % caxis([0 max(max(meanC))])
% set(h, 'EdgeColor', 'none');
% xlabel('Time (h)')
% ylabel('Depth (m)')
% cb.Label.String = 'Concentration (particles.m⁻¹)';
% title('MP - Aggregation off')
% 
% 
% f5 = figure(5); clf,
% h = pcolor(flipud(AggConcInter'));
% cb = colorbar;
% % caxis([0 max(max(meanC))])
% set(h, 'EdgeColor', 'none');
% xlabel('Time (h)')
% ylabel('Depth (m)')
% cb.Label.String = 'Concentration (particles.m⁻¹)';
% title('Agg - Aggregation off')
% 
% zF = zFinal{end};
% partI= cast(max(1, fix(zF/dz)), 'uint8'); % int array, index of each particle's current mesh
% partWs = VitesseNguyen(rhof, partFin, partI);
% f6 = figure(6); clf,
% scatter(partWs,[partFin.L],20,[partFin.Fb],'filled')
% xlabel('Velocity Ws (m.s⁻¹)')
% ylabel('Size (m)')
% cb = colorbar;
% cb.Label.String = 'Biological matter fraction Fb';
% 
% f7 = figure(7); clf,
% h = pcolor(flipud(MeanContent')); 
% set(h, 'EdgeColor', 'none');
% 
% scatter(1:P,-zFinal{end},20,AggContent{end},'filled')
% %%
% figures = [f1 f2 f3 f4 f5 f6];
% for iF = 1:length(figures)
%     savefig(figures(iF), [saveTo runID{:} '_f' num2str(iF) '.fig'])
%     exportgraphics(figures(iF), [saveTo runID{:} '_f' num2str(iF) '.png'])
% end


% f2 = figure(2); clf,
% p = zFinal{end};
% plot(-sort(p(nMP+1:end)),'.')
% hold on 
% plot(-sort(p(1:nMP)),'.')
% hold off
% 
% savefig(f2, [saveTo runID{:} '_f2.fig'])
% exportgraphics(f2, [saveTo runID{:} '_f2.png'])
% 
% %  while isempty(zFinal{end})
% %     zFinal = zFinal(1:end-1);
% %  end
%     %%
% pp = zFinal{end};
% ppSortMP = sort(pp(1:nMP));
% 
% topBound = dz;
% j = 1;
% count = zeros(N,1);
% 
% for i=1:length(ppSortMP)
%     while ppSortMP(i) > topBound
%         topBound = topBound+dz;
%         j = j+1;
%     end
%     count(j) = count(j)+1;
% end, clear i,
% CfMP = count/dz;
% clear ppSortMP topBound j,
% 
% ppSortAgg = sort(pp(nMP+1:end));
% 
% topBound = dz;
% j = 1;
% count = zeros(N,1);
% 
% for i=1:length(ppSortAgg)
%     while ppSortAgg(i) > topBound
%         topBound = topBound+dz;
%         j = j+1;
%     end
%     count(j) = count(j)+1;
% end, clear i,
% CfAgg = count/dz;
% clear ppSortAgg topBound j,
% 
% 
% f3 = figure(3); clf,
% plot(CfMP,-z_, CfAgg,-z_)
% 
% savefig(f3, [saveTo runID{:} '_f3.fig'])
% exportgraphics(f3, [saveTo runID{:} '_f3.png'])

end