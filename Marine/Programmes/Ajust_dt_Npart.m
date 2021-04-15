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
% nPart_test = [10000 20000 30000];
nPart = 50000;

tf = 1500;
%tf_test = [60 60*60 60*60*12 60*60*24];

% cKs = 1;
cKs_test = [0 0.1 0.5 1 5 10];
%cKs_test = [0];

results = {};
i = 0;
%for tf = tf_test
%for dt=dt_test
% for nPart = nPart_test
for cKs = cKs_test
    i = i+1;
    disp(class(nPart))
    disp(["nPart = ", num2str(nPart)])
    figure(1)
    [ti, MSEi] = Ajust_model(dt, nPart, tf, cKs);
    results(i,:) = {ti, MSEi};
end

leg = {};
inter = [];

f = figure(2); clf;
hold on
xlabel('Time (s)');
ylabel('Delta Concentration');
for i=1:size(results,1)
    
    if results{i,1}(end) == 1 && results{i,2}(end) == 1
        results{i,1} = results{i,1}(1:end-1);
        results{i,2} = results{i,2}(1:end-1);
    end
    
    plot(results{i,1}, sqrt(results{i,2}))
    xlim([0 tf])
    %leg(i) = {['tf = ',num2str(tf_test(i))]};
    %inter = [inter,'-', num2str(tf_test(i))];
    %leg(i) = {['dt = ',num2str(dt_test(i))]};
    %inter = [inter,'-', num2str(dt_test(i))];
%     leg(i) = {['nPart = ',num2str(nPart_test(i))]};
%     inter = [inter,'-', num2str(nPart_test(i))];
    leg(i) = {['cKs = ',num2str(cKs_test(i))]};
    inter = [inter,'-', num2str(cKs_test(i))];
end
legend(leg, 'Location','best');
hold off

% f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      inter,'_cKs', num2str(cKs),'.eps'];
% f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',inter,'_tf',...
%      num2str(tf),'_cKs', num2str(cKs),'.eps'];
% f_name = ['../../Ian/Results/nPart', inter,'_dt',num2str(dt),'_tf',num2str(tf),...
%     '_cKs', num2str(cKs),'.eps'];
f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
    num2str(tf),'_cKs', inter,'_var.eps'];
exportgraphics(f,f_name,'ContentType','vector');


% for i=1:size(results,1)
%     fi = figure(i+2);
%     clf
%     plot(results{i,1}, sqrt(results{i,2}))
%     xlim([0 tf])
%     legend(['nPart = ',num2str(nPart_test(i))],'Location','southeast');
%     xlabel('Time (s)');
%     ylabel('Delta Concentration (MPs.m⁻¹)');
%     fi_name = ['../../Ian/Results/nPart', num2str(mnPart),'-',num2str(nnPart),...
%         '-',num2str(MnPart),'_dt',num2str(dt),'_tf',num2str(tf),'_nPart',...
%         num2str(nPart_test(i)),'_cKs', num2str(cKs),'.eps'];
%     exportgraphics(fi,fi_name,'ContentType','vector');
% end


%  fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      inter,'_cKs', num2str(cKs),'.mat'];
% fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',inter,'_tf',...
%      num2str(tf),'_cKs', num2str(cKs),'.mat'];
% fvar_name = ['../../Ian/Results/nPart', inter,'_dt',num2str(dt),'_tf',...
%      num2str(tf),'_cKs', num2str(cKs),'.mat'];
fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
   num2str(tf),'_cKs', inter,'_var.mat'];
save(fvar_name, 'results')