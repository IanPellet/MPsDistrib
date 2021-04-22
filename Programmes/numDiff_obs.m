function [results] = numDiff_obs(VarObs, Ks, N_test)
%NUMDIFF_OBS Plots the error for different number of meshes
%
% VarObs is a string that can take 2 values : 'Ks' or 'Ws' corresponding to
% the variable parammeter
%
% Ks value(s) of diffusive turbidity (m².s⁻¹)
%
% Ws value(s) of speed (m.s⁻¹)
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
% elseif strcmp(VarObs, 'Ws')
%     Ws_test = Ws;
else
    disp('Invalid argument VarObs')
end
%path = '../../Ian/Results/'; 
path = './';
Ws = 0;
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
            erroriN = Transport_Eulerian(N, Ks);
            N_error(iN,1) = N;
            N_error(iN,2) = erroriN;
        end

        Emin = min(N_error(:,2)); % minimal error
        %imin = find(N_error(:,2) == Emin); % index of the Emin
        Nmin = N_error(N_error(:,2) == Emin); % N with minimal error
        %results(iKs,:) = {Ks, Nmin, Emin, N_error}; % store results
        results.bilan(iKs,:) = {Ks, Nmin, Emin}; % store results
        results.Ks(iKs) = Ks; % store results
        results.N = N_error(:,1); % store results
        results.Error(iKs,:) =  N_error(:,2); % store results

        Ks_Nmin(iKs,1) = Ks;
        Ks_Nmin(iKs,2) = Nmin;

        inter = [inter '-' num2str(Ks)];
        
        %% Plot error
        fError = figure(1); clf;
        semilogy(N_error(:,1),N_error(:,2))
        xlabel('N');
        ylabel('Mean Square Error');
        ttl = ['Ks = ' num2str(Ks) 'm².s⁻¹ ; Ws = ' num2str(Ws) 'm.s⁻¹'];
        title(ttl)
        
        fname = ['Ks', num2str(Ks),'_Ws', num2str(Ws),'_N',...
            num2str(min(N_test)),'-',num2str(length(N_test)),'-',...
            num2str(max(N_test))];
%         exportgraphics(fError,[path, 'ErrorKs_', fname,'.eps'],'ContentType','vector');
%         savefig(fError,[path, 'ErrorKs_', fname,'.fig']);
    end   
end

%% Plot errors on the same fig
% fError2 = figure(3); clf;
% for j = 1:length(Ks_test)
%     j_name (j)= ['Ks = ' num2str(Ks_test(j)) 'm².s⁻¹'];
% end
%     disp(j_name)
%     semilogy(results.N,results.Error,'DisplayName',j_name)
semilogy(results.N,results.Error)
hold on
leg = {};
for i = 1:length(results.Ks)
    li = ['Ks = ' num2str(Ks_test(i)) 'm².s⁻¹'];
    leg{i} = li;
end
legend(leg, 'Location', 'best');
xlabel('N');
ylabel('Mean Square Error');
hold off

%% Plot min results
fNmin = figure(2); clf;

if strcmp(VarObs, 'Ks')
    pcolor(results.Ks,results.N,results.Error');
    hold on
    plot(Ks_Nmin(:,1), Ks_Nmin(:,2), 'm')
    xlabel('Ks (m².s⁻¹)');
    ylabel('N_m_i_n');
    hold off
end


%% Save results

if strcmp(VarObs, 'Ks')
    
    files = ['Ks', num2str(min(Ks_test)),'-',...
        num2str(length(Ks_test)),'-', num2str(max(Ks_test)),'_Ws',...
        num2str(Ws),'_N',num2str(min(N_test)),'-',num2str(length(N_test)),...
        '-', num2str(max(N_test))];
   
end

save([path,files,'.mat'], 'results', [VarObs '_Nmin'])

%% Save plot
% exportgraphics(fError2,[path, 'ErrorUnit_', files, '.eps'],'ContentType','vector');
% savefig(fError2,[path, 'ErrorUnit_', files, '.fig']);
% exportgraphics(fNmin,[path, 'Nmin_', files, '.eps'],'ContentType','vector');

