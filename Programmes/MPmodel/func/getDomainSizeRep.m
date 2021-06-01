function [meanSize, stdSize, s_, s] = getDomainSizeRep(topBound, bottomBound, mp, zPart)

iDomain = cell2mat(zPart) >= topBound & cell2mat(zPart) <= bottomBound;

sizeHistory = cell(size(zPart));
for i=1:length(zPart)
    sizeHistory{i} = [mp(iDomain(i,:)).size_];
end, clear i,

ds = 5e-6;
s = 0:ds:600e-6;
s_ = s(1:end-1)+ds/2;
[meanSize, stdSize] = getMeanConc(sizeHistory, length(s_), ds);

% figure
% hold on
% plot(s_, meanSize, 'b')
% plot(s_, meanSize+2*stdSize, 'b--')
% plot(s_, meanSize-2*stdSize, 'b--')
% hold off
end