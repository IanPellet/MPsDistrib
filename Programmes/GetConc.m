function [Conc] = GetConc(zHistory, iStart, iEnd,N, dz)
% GETCONC Get particle concentration from particle position
%
% zHistory : particle position during simulation
% iStart : 
% iEnd : 
% N : 
% dz : water column discretisation step (m)

    Conc = NaN(length(zHistory),N);
    for hStep = 1:length(zHistory)

        pp = zHistory{hStep};
        ppSort = sort(pp(iStart:iEnd));

        topBound = dz;
        j = 1;
        count = zeros(N,1);

        for i=1:length(ppSort)
            while ppSort(i) > topBound
                topBound = topBound+dz;
                j = j+1;
            end
            count(j) = count(j)+1;
        end, clear i,
        Conc(hStep,:) = count/dz;
    end
end