function [results] = numDiff_obs(VarObs, Ks, Ws, N_test)
%NUMDIFF_OBS Plots the error for different number of meshes
%
% VarObs is a string that can take 2 values : 'Ks' or 'Ws' corresponding to
% the variable parammeter
%
% Ks value(s) of diffusive turbidity (m².s⁻¹)
%
% Ws value(s) of spedd (m.s⁻¹)
%
% N_test values taken for the number of meshes
%
%

%% Variables init
% N_test = [10 50 100 200 500]; 
%N_tested = 1000:100:2000;
% N_test = 10:10:500;
% N_test = linspace(10,510,11);
if strcmp(VarObs, 'Ks')
    Ks_test = Ks;
elseif strcmp(VarObs, 'Ws')
    Ws_test = Ws;
else
    disp('Invalid argument VarObs')
end
path = ['../../Ian/Results/'];
%v_test = [1e-4 1e-3 1e-2 1e-1];
%v_test = linspace(1e-2, 1e-1, 4);
%v = 1e-3;

%% Data structures init
results = {};
inter = [];



%% Run simulations
if strcmp(VarObs, 'Ks')
    Ks_Nmin = zeros(length(Ks_test),2);
    iKs = 0;
    for Ks = Ks_test
        iKs = iKs+1;

        N_error = zeros(length(N_test),2);
        iN = 0;
        for N = N_test
            iN = iN+1;
            erroriN = Transport_Eulerian(N, Ks, Ws);
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
        
        %% Plot error
        fError = figure(1); clf;
        semilogy(N_error(:,1),N_error(:,2))
        xlabel('N');
        ylabel('Error (mps.m⁻³)');
        ttl = ['Ks = ' num2str(Ks) 'm².s⁻¹ ; Ws = ' num2str(Ws) 'm.s⁻¹'];
        title(ttl)
        
        fname = ['Ks', num2str(Ks),'_Ws', num2str(Ws),'_N',...
            num2str(min(N_test)),'-',num2str(length(N_test)),'-',...
            num2str(max(N_test))];
        exportgraphics(fError,[path, 'ErrorKs_', fname,'.eps'],'ContentType','vector');
    end

    
elseif strcmp(VarObs, 'Ws')
    Ws_Nmin = zeros(length(Ws_test),2);
    iWs = 0;
    for Ws = Ws_test
        iWs = iWs+1;

        N_error = zeros(length(N_test),2);
        iN = 0;
        for N = N_test
            iN = iN+1;
            erroriN = Transport_Eulerian(N, Ks, Ws);
            N_error(iN,1) = N;
            N_error(iN,2) = erroriN;
        end

        Emin = min(N_error(:,2)); % minimal error
        imin = find(N_error(:,2) == Emin); % index of the Emin
        Nmin = N_error(imin); % N with minimal error
        results(iWs,:) = {Ws, Nmin, Emin, N_error}; % store results

        Ws_Nmin(iWs,1) = Ws;
        Ws_Nmin(iWs,2) = Nmin;

        inter = [inter '-' num2str(Ws)];
        
        %% Plot error
        fError = figure(1); clf;
        semilogy(N_error(:,1),N_error(:,2))
        xlabel('N');
        ylabel('Error (mps.m⁻³)');
        ttl = ['Ks = ' num2str(Ks) 'm².s⁻¹ ; Ws = ' num2str(Ws) 'm.s⁻¹'];
        title(ttl)
        
        fname = ['Ks', num2str(Ks),'_Ws', num2str(Ws),'_N',...
            num2str(min(N_test)),'-',num2str(length(N_test)),'-',...
            num2str(max(N_test))];
        exportgraphics(fError,[path,'ErrorWs_',fname,'.eps'],'ContentType','vector');
    end    
end

%% Plot errors on the same fig
fError2 = figure(3); clf;
for j = 1:length(Ks_test)
    j_error = results{j,4};
    j_name = ['Ks = ' num2str(Ks_test(j)) 'm².s⁻¹'];
    disp(j_name)
    semilogy(j_error(:,1),j_error(:,2),'DisplayName',j_name)
    hold on
end
legend('Location','best');
xlabel('N');
ylabel('Error (mps.m⁻³)');
hold off

%% Plot min results
fNmin = figure(2); clf;

if strcmp(VarObs, 'Ks')
    plot(Ks_Nmin(:,1), Ks_Nmin(:,2))
    xlabel('Ks (m².s⁻¹)');
    ylabel('N_m_i_n');
        
elseif strcmp(VarObs, 'Ws')
    plot(Ws_Nmin(:,1), Ws_Nmin(:,2))
    xlabel('Ws (m.s⁻¹)');
    ylabel('N_m_i_n');
end


%% Save results

if strcmp(VarObs, 'Ks')
    files = ['Ks', num2str(min(Ks_test)),'-',...
        num2str(length(Ks_test)),'-', num2str(max(Ks_test)),'_Ws',...
        num2str(Ws),'_N',num2str(min(N_test)),'-',num2str(length(N_test)),...
        '-', num2str(max(N_test))];
        
elseif strcmp(VarObs, 'Ws')
    files = ['Ws', inter,'_Ks', num2str(Ks),'_N',...
    num2str(min(N_test)),'-',num2str(length(N_test)),'-', num2str(max(N_test))];
end

save([files,'.mat'], 'results', [VarObs '_Nmin'])

%% Save plot
exportgraphics(fError2,[path, 'Error2_', files, '.eps'],'ContentType','vector');
exportgraphics(fNmin,[path, 'Nmin_', files, '.eps'],'ContentType','vector');

