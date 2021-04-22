dotMatFiles = struct2cell(dir(fullfile(['../../Ian/Results/MSE_alpha/0412_mse-a_0.5-2.5_10_mina_*.mat'])))';

N = size(dotMatFiles,1);
min_a_array = zeros(N);
for i = 1:N
    file = fullfile(['../../Ian/Results/MSE_alpha/',dotMatFiles{i}]);
    storedStructure = load(file, 'min_a_Ks'); % Load in ONLY the myVar variable.
    min_a_array(i,:) = storedStructure.min_a_Ks;  % Assign it to a new variable with different name.
    clear('storedStructure'); % If it's really not needed any longer.
end