function [avgCDF, stdCDF, avgNDom, s] = getDomainSizeCDF(topBound, bottomBound, mp, zPart)
%GETDOMAINSIZECDF Get CDF of particles sizes in the domain def by bounds
%
% INPUTS :
% topBound : double, domain's top boundary (m)
% bottomBound : double, domain's bottom boundary (m)
% mp : MP object array, list of MPs in the water column
% zPart : cell array of double array, position of MPs at each time step (m)
%
% OUTPUTS :
% meanCDF : double array, time average on the time steps in zPart of the
% particle size CDF (m)
% stdCDF : double array, time standard deviation on the time steps in zPart
% of the particle size CDF (m)
% avgNDom : int, time average number of MPs in the domain on the time steps
% in zPart
% s : double array, particle sizes at which the CDF is computed (m)
%

% Get indexes of particles in the domain
iDomain = cell2mat(zPart) >= topBound & cell2mat(zPart) <= bottomBound;
nPartDom = sum(iDomain,2); % number of particles in the domain at each time step
avgNDom = mean(nPartDom); % time avg of number of part in the domain

% Get the sizes of the particles in the domain at each time step
sizeHistory = cell(size(zPart));
for i=1:length(zPart)
    sizeHistory{i} = [mp(iDomain(i,:)).size_];
end, clear i,

% Define the size interval to evaluate CDF
ds = 5e-6;
s = 0:ds:max([mp.size_]);

% Get size CDF at each time step
cdfHistory = cell(size(sizeHistory));
for hStep = 1:size(sizeHistory,1)
    pp = sizeHistory{hStep};
    ppSort = sort(pp);
    if length(ppSort) > 1
        x = [1:length(ppSort)]/length(ppSort);
        cdfHistory{hStep} = interp1(ppSort,x,s);
    end
end, clear hStep pp ppSort x,

cdfMat = cell2mat(cdfHistory);

% Compute average and std of the size CDF computed
avgCDF = mean(cdfMat);
stdCDF = std(cdfMat);

end