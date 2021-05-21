function [mpFinal, zFinal] = MP_fragmentation(mpInit, zInit, iMP)
%MP_FRAGMENTATION Creates a new MP array and position array with fragments
%of MP at indexes iMP
%
% PARAMETERS
% mpInit : MP objects array, modeled particles
% zInit : double array, particle's position (m)
% iMP : logical array, MP to fragment in mpInit
% 
% OUTPUTS
% mpFinal : MP objects array, modeled particles with created fragments
% zFinal : double array, particle's position with created fragments (m)
%

%% Create Fragments
frag(2,sum(iMP)) = MP;
% Get particle to fragment
mpInit0 = mpInit(iMP);
for imp = 1:length(mpInit0)
    mp0 = mpInit0(imp);
    % Create two fragments 
    frag(2,:) = MP(mp0.size_/2, mp0.rho_, mp0.rhow_, 0);
    frag(1,:) = MP(mp0.size_/2, mp0.rho_, mp0.rhow_, 0);
end, clear imp mp0 mpInit0;

%% Fill returned arrays
mpInter = mpInit;
mpInter(iMP) = frag(1,:);
mpFinal = [mpInter frag(2,:)];

zFinal = [zInit zInit(iMP)];

end