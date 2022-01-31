% Ploting size distribution at the top of the column for different data
Tcozart = readtable('../Data/cozartSizeDist.txt');
load('../Data/Emilie/allData.mat')
Temilie = Tdata; clear Tdata,
Temilie = Temilie(Temilie(:,'depth').Variables<=1,:);
Tanaelle = readtable('../Data/data_mps.txt');
Tanaelle = Tanaelle(Tanaelle(:,'depth').Variables<=1,:);
Tanaelle = Tanaelle(Tanaelle(:,'type').Variables~=0,:);

bins = [Tcozart(:,'m').Variables ; Tcozart(end,'M').Variables];
sumSize = cumsum(Tcozart(:,'N').Variables);
NsizeNorm = Tcozart(:,'N').Variables./Tcozart(:,'dS').Variables;

Sizes = Tcozart(:,'nomSize').Variables;
SizesAnaelle = sort(unique(Tanaelle(:,'med_size').Variables))*1e-3;

BinsAnaelle = [50,200,250,500,1000:500:5000,30000]*1e-3;
dSAnaelle = diff(BinsAnaelle);

derivS = diff(log(sumSize))./diff(log(Sizes));

SizeArea = Temilie(Temilie(:,'type').Variables == 'Area','size').Variables*1e3;
% SizeArea = Temilie(:,'size').Variables*1e3;

countSizes = histcounts(SizeArea,bins);
NsN = countSizes./Tcozart(:,'dS').Variables';

cSizesAnaelle = zeros(size(SizesAnaelle));
for iS = 1:length(SizesAnaelle)
    cSizesAnaelle(iS) = sum(Tanaelle(Tanaelle(:,'med_size').Variables==SizesAnaelle(iS)*1e3,'n').Variables);
end, clear iS,
NsNAnaelle = cSizesAnaelle./dSAnaelle';


% Fit
% cfCozart = fit(Sizes(15:end),NsizeNorm(15:end),'Power1');
% cfEmilie = fit(Sizes(NsN~=0),NsN(NsN~=0)','Power1');
% cfAnaelle = fit(SizesAnaelle(NsNAnaelle~=0)*1e-3,NsNAnaelle(NsNAnaelle~=0),'Power1');

toFit = {[Sizes(14:end) NsizeNorm(14:end)] ; [Sizes(NsN~=0) NsN(NsN~=0)'] ; [SizesAnaelle(6:end) NsNAnaelle(6:end)]};
FittedCoef = cell(size(toFit));
for iF = 1:length(toFit)
    coefPlot = toFit{iF};
    xfit = log(coefPlot(:,1));
    yfit = log(coefPlot(:,2));
    X = [ones(length(xfit),1) xfit];
    FittedCoef{iF} = X\yfit;
end, clear iF,


figure(1), clf,
loglog(Sizes,NsizeNorm,'+','DisplayName','Cozart')
hold on
loglog(Sizes,NsN,'+','DisplayName','Emilie')
loglog(SizesAnaelle,NsNAnaelle,'+','DisplayName','Anaelle')

DispFit = {'Fit Cozart' ; 'Fit Emilie' ; 'Fit Anaelle'};
for iF = 1:length(FittedCoef)
    coefPlot = FittedCoef{iF};
    yy = exp(coefPlot(1))*Sizes.^coefPlot(2);
    loglog(Sizes,yy,'DisplayName',DispFit{iF})
end, clear iF,

hold off
legend
xlabel('Size (mm)')
ylabel('Normalized abundance')


cumSumData = {flipud(cumsum(flipud(Tcozart(:,'N').Variables))) ;
    flipud(cumsum(flipud(countSizes'))) ;
    flipud(cumsum(flipud(cSizesAnaelle)))};
sizeData = {Sizes ; Sizes ; SizesAnaelle};
fitDataI = {find(Sizes>=5, 1 ):length(sizeData{1})-1 ; cumSumData{2}~=0 ; 1:length(sizeData{3})-1};
FittedCoefCS = cell(size(cumSumData));
for iF = 1:length(cumSumData)
    sizeFit = sizeData{iF};
    NFit = cumSumData{iF};
    xfit = log(sizeFit(fitDataI{iF}));
    yfit = log(NFit(fitDataI{iF}));
    X = [ones(length(xfit),1) xfit];
    FittedCoefCS{iF} = X\yfit;
end, clear iF,


figure(2), clf,
DispName = {'data_C' ; 'data_S' ; 'data_A'};
DispFit = {'Fitted D1' ; 'Fitted D2' ; 'Fited D3'};
coolors = {'#540D6E' ; '#0B5D1E' ; '#00A6A6'};
for iF = 1:length(FittedCoefCS)
    loglog(sizeData{iF},cumSumData{iF},'x','DisplayName',DispName{iF},'color', coolors{iF}, 'MarkerSize', 7)
    hold on
    coefPlot = FittedCoefCS{iF};
    yy = exp(coefPlot(1))*[SizesAnaelle(1) ; Sizes].^coefPlot(2);
    loglog([SizesAnaelle(1) ; Sizes],yy,'DisplayName',...
        [num2str(exp(coefPlot(1))) '*L^-^D ; D = ' num2str(-coefPlot(2))] ,'color', coolors{iF},'lineWidth',2)
end, clear iF,

hold off
legend('Location','Southwest')
xlabel('L (mm)')
ylabel('Cumulated abundance')

