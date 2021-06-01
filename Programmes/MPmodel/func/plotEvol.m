saveInt = [86400 172800 259200 345600 432000];
filePath = '/media/ian/Transcend/MPsDistrib/Results/MP_runStabTest/Data/2021-6-1-9-5-28_691382-zHist';
dtStab = 1800; 
dt = 10;
tf = 432000;
dCinter = 0:dtStab:tf;
dCinter_ = dCinter(1:end-1)+dtStab/2;
timeLine = 0:dt:tf;
nCase = dtStab/dt;
L = 67.3456;
N = 50;
dz = L/N;
z=0:dz:L; % z : boundaries of the meshes
z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh  

i = 0;
C = cell(length(dCinter)-1,2);
for iload = 1:length(saveInt)
    file = load([filePath num2str(saveInt(iload)), '.mat'], "zHistory");
%     history = [file.zHistory];
    history = zHistory;
    if isempty(history{end})
        history = history(1:end-1);
    end
    iStart = 0;
    iEnd = 0;

    for i30 = 1:length(history)/nCase
        iStart = iEnd + 1;
        iEnd = iEnd + nCase;
        i = i+1;
        [C{i,1}, C{i,2}] = getMeanConc(history(iStart:iEnd), N, dz);
    end
    clear hist i30,
end, clear i iload file history,



meanC = cell2mat(C(:,1));
pcolor(dCinter_',-z_,meanC')
colorbar