load('../Results/sizeRep/DataSim.mat')
ID = strrep(join(split(num2str(clock)),'-'),'.','_'); % unique ID to name files
ID = ID{:};

sizes = [2e-6 1e-3]; % Particle sizes

% Size array init
mpSize = ones(nPart,1);
iMP = 0;
nMPsize = round(nPart/length(sizes));
for iSize = 1:length(sizes)
    iMP = iMP+1;
    iMPend = iMP + nMPsize;
    if iMPend > nPart
        iMPend = nPart;
    end
    mpSize(iMP:iMPend) = mpSize(iMP:iMPend)*sizes(iSize); 
end, clear iSize,

% MP init
mp = getMPlist(nPart, mpSize, rhop, rhow, frag);
clear mpSize,

% Initial position
zInit = ones(nPart,1)*L/2;

[zFinal,dt] = MP_simulator(mp, zInit, K, dK, L, dz, tf, dt_test, tavg);
[meanConc, stdConc] = getMeanConc(zFinal, length(z_), dz);

plot(meanConc, -z_);


% save(['../Results/sizeRep/DataSim' ID '.mat'], 'ID', 'sizes', 'mp', 'zInit',...
%     'zFinal', 'dt', 'meanConc', 'stdConc')