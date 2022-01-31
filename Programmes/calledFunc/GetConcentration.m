function [hConc] = GetConcentration(zPart, nCat, dz)
%% GETCONCENTRATION  Computes particle concentration from there positions
% 
% PARAMETERS
% zPart : double array, particle positions
% nCat : int, number of meshes to discretize the water column
% dz : double, space discretization step (m)

    hConc = NaN(size(zPart,1),nCat);
    for hStep = 1:size(zPart,1)

        pp = zPart{hStep};

        ppSort = sort(pp);

        topBound = dz;
        j = 1;
        count = zeros(nCat,1);

        for i=1:length(ppSort)
            if ppSort(i) > topBound
                topBound = topBound+dz;
                j = j+1;
            end
            count(j) = count(j)+1;
        end, clear i,

        hConc(hStep,:) = count'/dz;

    end, clear hStep pp histi,
end