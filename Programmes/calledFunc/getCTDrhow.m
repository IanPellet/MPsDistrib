function [rhow] = getCTDrhow(day,z)
%GETCTDRHOW Get water density from CTD data
%
% PARAMETERS
% day : string, sampling day {'10fev', '3fev'}
% z : double array, water column discrtisation (m)
%
% OUTPUT
% rhow : double array, water density on water column discrtisation (kg.m⁻³)
%
    
    if strcmp(day,'10fev')
        CTDfile= '../Data/CTD_Marine/CTD7708-20200210-0841-ave'; 
    elseif strcmp(day,'3fev')
        CTDfile = '../Data/CTD_Marine/CTD7708-20200203-0849';
    else
        error('Unrecognised day');
    end
    
    CTD = load(CTDfile); % load file
    CTD = sortrows(CTD); % sort data along depth
    [~,i,~] = unique(CTD(:,1)); % filter unique depth data
    dens = CalculDensite(CTD(i,3),CTD(i,2)); % compute density
    rhow = interp1(CTD(i,1),dens,z,'pchip'); % interpolation on water column discrtisation 
end