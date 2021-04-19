%NUMDIFF_OBS Plots the error for different number of meshes

%% Variables init
N_test = linspace(10,510,11);

Ks_test = [0.001 0.01 0.1 1];

%v_test = [1e-4 1e-3 1e-2 1e-1];
%v_test = linspace(1e-2, 1e-1, 4);
v = 1e-3;

%% Data structures init
results = {};
Ks_Nmin = zeros(length(Ks_test),2);
inter = [];

%% Plot error init
fError = figure(1); clf;
hold on
xlabel('N');
ylabel('DeltaC (%)');

%% Run simulations
iKs = 0;
for Ks = Ks_test
    iKs = iKs+1;
    
    N_error = zeros(length(N_test),2);
    iN = 0;
    for N = N_test
        iN = iN+1;
        erroriN = Transport_Eulerian(N, Ks, v);
        N_error(iN,1) = N;
        N_error(iN,2) = erroriN;
    end
    
    Emin = min(N_error(:,2)); % minimal error
    imin = find(N_error(:,2) == Emin); % index of the Emin
    Nmin = N_error(imin); % N with minimal error
    results(iKs,:) = {Ks, Nmin, Emin, N_error}; % store results
    
    Ks_Nmin(iKs,1) = Ks;
    Ks_Nmin(iKs,2) = Nmin;
    
    inter = [inter '-' num2str(Ks)];
    
    plot(N_error(:,1),N_error(:,2)*100, 'DisplayName',['Ks = ' num2str(Ks) 'm².s⁻¹'])
end
legend('Location','best');
hold off


%% Plot min results
fNmin = figure(2); clf;

plot(Ks_Nmin(:,1), Ks_Nmin(:,2))

xlabel('Ks (m².s⁻¹)');
ylabel('N_m_i_n');

%% Save results
files = ['../../Ian/Results/Ks', inter,'_v', num2str(v),'_N',...
    num2str(min(N_test)),'-',num2str(length(N_test)),'-', num2str(max(N_test))];

save([files,'.mat'], 'results', 'Ks_Nmin')

%% Save plots
exportgraphics(fError,[files,'_Error.eps'],'ContentType','vector');
exportgraphics(fNmin,[files,'_Nmin.eps'],'ContentType','vector');

