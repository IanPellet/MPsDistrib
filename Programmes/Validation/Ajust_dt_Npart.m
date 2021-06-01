
% mdt = 0.01;
% Mdt = 5;
% ndt = 10;
% dt_test = linspace(mdt, Mdt, ndt);
% dt_test = [1 5 10 50 100];
% dt = 1;

% mnPart = 50;
% MnPart = 10000;
% nnPart = 6;
% nPart_test = linspace(mnPart, MnPart, nnPart);
% nPart_test = [50 100 500 1000 5000 10000 50000];
nPart = 50e3;

tf = 1e5;
% tf_test = [60*60*24];

% Ks = 0.01;
Ks_test = [0.001 0.01 0.1 1];
%Ks_test = [0];

%v_test = [1e-4 1e-3 1e-2 1e-1];
%v_test = linspace(1e-2, 1e-1, 4);
v = 1e-3;

dtest = 60*60*2;

%dt = fix(1/max(v_test));
dt = fix(1/max(v));
disp(["dt = ",num2str(dt),"s"])

results = {};
i = 0;
% for tf = tf_test
% for dt=dt_test
% for nPart = nPart_test
for Ks = Ks_test
% for v = v_test
    i = i+1;
    disp(class(nPart))
    disp(["nPart = ", num2str(nPart)])
    %figure(1)
    [ti, errori] = Ajust_model(dt, nPart, tf, Ks, v, dtest);
    results(i,:) = {ti, errori};
end

leg = {};
inter = [];

f = figure(2); clf;
hold on
xlabel('Time (s)');
ylabel('DeltaC (%)');
for i=1:size(results,1)
    
    if results{i,1}(end) == 1 || results{i,2}(end) == 1
        results{i,1} = results{i,1}(1:end-1);
        results{i,2} = results{i,2}(1:end-1);
    end
    
    plot(results{i,1}, results{i,2}*100)
    xlim([0 tf])
%     leg(i) = {['tf = ',num2str(tf_test(i))]};
%     inter = [inter,'-', num2str(tf_test(i))];
%     leg(i) = {['dt = ',num2str(dt_test(i)),'s']};
%     inter = [inter,'-', num2str(dt_test(i))];
%     leg(i) = {['nPart = ',num2str(nPart_test(i))]};
%     inter = [inter,'-', num2str(nPart_test(i))];
    leg(i) = {['Ks = ',num2str(Ks_test(i)), 'm².s⁻¹']};
    inter = [inter,'-', num2str(Ks_test(i))];
%     leg(i) = {['v = ',num2str(v_test(i)),'m.s⁻¹']};
%     inter = [inter,'-', num2str(v_test(i))];

end
legend(leg, 'Location','best');
hold off

% f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      inter,'_Ks', num2str(Ks),'_v', num2str(v),'_dtest', num2str(dtest),'_error.eps'];
% f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',inter,'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_v', num2str(v),'_dtest', num2str(dtest),'_error.eps'];
% f_name = ['../../Ian/Results/nPart', inter,'_dt',num2str(dt),'_tf',num2str(tf),...
%     '_Ks', num2str(Ks),'_v', num2str(v),'_dtest', num2str(dtest),'_error.eps'];
f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
    num2str(tf),'_Ks', inter,'_v', num2str(v),'_dtest', num2str(dtest),'_error.eps'];
% f_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_v', inter,'_dtest', num2str(dtest),'_error.eps'];

exportgraphics(f,f_name,'ContentType','vector');



% fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      inter,'_Ks', num2str(Ks),'_v', num2str(v),'_dtest', num2str(dtest),'_error.mat'];
% fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',inter,'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_v', num2str(v),'_dtest', num2str(dtest),'_error.mat'];
% fvar_name = ['../../Ian/Results/nPart', inter,'_dt',num2str(dt),'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_v', num2str(v),'_dtest', num2str(dtest),'_error.mat'];
fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
   num2str(tf),'_Ks', inter,'_v', num2str(v),'_dtest', num2str(dtest),'_error.mat'];
% fvar_name = ['../../Ian/Results/nPart', num2str(nPart),'_dt',num2str(dt),'_tf',...
%      num2str(tf),'_Ks', num2str(Ks),'_v', inter,'_dtest', num2str(dtest),'_error.mat'];
save(fvar_name, 'results')