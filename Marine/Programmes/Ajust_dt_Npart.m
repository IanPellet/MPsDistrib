%mdt = 1;
%Mdt = 5;
%ndt = 10;
%dt_test = linspace(mdt, Mdt, ndt);
dt = 1;

mnPart = 50;
MnPart = 10000;
nnPart = 6;
%nPart_test = linspace(mnPart, MnPart, nnPart);
nPart_test = [50 100 500 1000 5000 10000];

tf = 1000;

results = {};
i = 0;
%for dt=dt_test
for nPart = nPart_test
    i = i+1;
    disp(class(nPart))
    disp(["nPart = ", num2str(nPart)])
    figure(1)
    [ti, MSEi] = Ajust_model(dt, nPart, tf);
    results(i,:) = {ti, MSEi};
end

leg = {};

f = figure(2); clf;
hold on
xlabel('Time (s)');
ylabel('Delta Concentration (MPs.m⁻¹)');
for i=1:size(results,1)
    
    if results{i,1}(end) == 1 && results{i,2}(end) == 1
        results{i,1} = results{i,1}(1:end-1);
        results{i,2} = results{i,2}(1:end-1);
    end
    
    plot(results{i,1}, sqrt(results{i,2}))
    xlim([0 tf])
    leg(i) = {num2str(nPart_test(i))};
end
legend(leg, 'Location','southeast');
hold off

f_name = ['../../Ian/Results/nPart', num2str(mnPart),'-',num2str(nnPart),'-',num2str(MnPart),'_dt',num2str(dt),'_tf',num2str(tf),'.eps'];
exportgraphics(f,f_name,'ContentType','vector');


for i=1:size(results,1)
    fi = figure(i+2);
    clf
    plot(results{i,1}, sqrt(results{i,2}))
    xlim([0 tf])
    legend(['nPart = ',num2str(nPart_test(i))],'Location','southeast');
    xlabel('Time (s)');
    ylabel('Delta Concentration (MPs.m⁻¹)');
    fi_name = ['../../Ian/Results/nPart', num2str(mnPart),'-',num2str(nnPart),'-',num2str(MnPart),'_dt',num2str(dt),'_tf',num2str(tf),'_nPart',num2str(nPart_test(i)),'.eps'];
    exportgraphics(fi,fi_name,'ContentType','vector');
end

fvar_name = ['../../Ian/Results/nPart', num2str(mnPart),'-',num2str(nnPart),'-',num2str(MnPart),'_dt',num2str(dt),'_tf',num2str(tf),'.mat'];
save(fvar_name, 'results')