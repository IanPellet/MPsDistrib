path = '../../Ian/Results/varMP/';

wind = 10; % km/h
month = 02;

D_test = [750 1250 1750 2250]*1e-6;
% rhop_test = linspace(850,1300, 10);
rhop_test = [1010.5]; % kg.m⁻³
nPart = 50e3;

dt = 10;
tf = 1e5;
dt_test = 60*60;

results = cell(length(D_test)*length(rhop_test),4); % Memory allocation for storing results

i = 0;
for D = D_test
    for rhop = rhop_test
        i = i+1;
        [Zi, Ci] = varMP_model(D, rhop, nPart, dt, tf, dt_test, wind, month, path);
        results(i,:) = {D, rhop, Zi, Ci}; % Store values
    end
end


file_name = ['results_D',num2str(min(D)),'-',num2str(length(D)),'-',num2str(max(D)), '_rhop',num2str(rhop),...
    '_nPart',num2str(nPart), '_dt', num2str(dt),...
    '_tf', num2str(tf), '_dtest', num2str(dt_test),'_wind', num2str(wind),...
    '_month', num2str(month)];
save([path,file_name,'.mat'],'results');
