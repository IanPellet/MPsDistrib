function [StabC,tStabC] = getStability(meanConc, testBetween, dtStab, tC)
%GETSTABILITY computes the root mean square error between avg concentration
% profiles separated by testBetween
%
% Parameters :
% meanConc : cell arrays of each time average concentration profile
% testBetween :time between two profiles compared (s)
% tC : double array, times at which we have a value meanConc (s)
% 
% Returns :
% StabC : double array, root mean square error of profiles separated by
% testBetween seconds
% tStabC : double array, times at which we have a value for StabC (s)

% number of indexes corresponding to testBetween in meanConc
diStab = testBetween/dtStab; 

% Times at which we compute a value for StabC
tStabC = (tC(1:end-diStab)+tC(1+diStab:end))/2;


StabC = NaN(1,length(meanConc)-diStab);
for idC = 1:length(meanConc)-diStab
    % Compute root mean square error
    StabC(idC) = sqrt(mean((meanConc{idC}-meanConc{idC+diStab}).^2));
end, clear idC,

end