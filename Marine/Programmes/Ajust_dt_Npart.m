% mdt = 0.01;
% Mdt = 5;
% ndt = 10;
% dt_test = linspace(mdt, Mdt, ndt);
% dt_test = [0.1 1 10];
dt = 1;

% mnPart = 50;
% MnPart = 10000;
% nnPart = 6;
% nPart_test = linspace(mnPart, MnPart, nnPart);
% nPart_test = [50 100 500 1000 5000 10000 50000];
nPart = 50e3;

tf = 1000;
% tf_test = [60*60*24];

% Ks = 0.01;
Ks_test = [0 0.001 0.005 0.01 0.05 0.1 0.5 1];
%Ks_test = [0];

results = {};
i = 0;
% for tf = tf_test
% for dt=dt_test
% for nPart = nPart_test
for Ks = Ks_test
    i = i+1;
    disp(class(nPart))
    disp(["nPart = ", num2str(nPart)])
    %figure(1)
    [ti, errori] = Ajust_model(dt, nPart, tf, Ks);
    results(i,:) = {ti, errori};
end

leg = {};
inter = [];

f = figure(2); clf;
hold on
xlabel('Time (s)');
ylabel('Difference with initial condition (%)');
for i=1:size(results,1)
    
    if results{i,1}(end) == 1 && results{i,2}(end) == 1
        results{i,1} = results{i,1}(1:end-1);
        results{i,2} = results{i,2}(1:end-1);
    end
    
    plot(results{i,1}, results{i,2}*100)
    xlim([0 tf])
%     leg(i) = {['tf = ',num2str(tf_test(i))]};
%     inter = [inter,'-', num2str(tf_test(i))];
%     leg(i) = {['dt = ',num2str(dt_test(i))]};
%     inter = [inter,'-', num2str(dt_test(i))];
%     leg(i) = {['nPart = ',num2str(nPart_test(i))]};
%     inter = [inter,'-', num2str(nPart_test(i))];
    leg(i) = {['Ks = ',num2str(Ks_test(i))]};
    inter = [inter,'-', num2str(Ks_test(i))];
end
legend(leg, 'Location','best');
hold off

% f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      inter,'_Ks', num2str(Ks),'_Kcst.eps'];
% f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',inter,'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_Kcst.eps'];
% f_name = ['../../Ian/Results/nPart', inter,'_dt',num2str(dt),'_tf',num2str(tf),...
%     '_Ks', num2str(Ks),'_Kcst.eps'];
f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
    num2str(tf),'_Ks', inter,'_Kcst.eps'];
exportgraphics(f,f_name,'ContentType','vector');



% fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      inter,'_Ks', num2str(Ks),'_Kcst.mat'];
% fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',inter,'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_Kcst.mat'];
% fvar_name = ['../../Ian/Results/nPart', inter,'_dt',num2str(dt),'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_Kcst.mat'];
fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
   num2str(tf),'_Ks', inter,'_Kcst.mat'];
save(fvar_name, 'results')