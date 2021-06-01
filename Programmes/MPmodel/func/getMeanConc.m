function [meanConc, stdConc] = getMeanConc(zPart, z, z_, dz, L)
% zPart = zFinal(end-10:end);
    hConc = NaN(size(zPart,1),length(z_));
    for hStep = 1:size(zPart,1)

        pp = zPart{hStep};

        ppSort = sort(pp);

        topBound = dz;
        j = 1;
        count = zeros(size(z_));

        for i=1:length(ppSort)
            if ppSort(i) > topBound
                topBound = topBound+dz;
                j = j+1;
            end
            count(j) = count(j)+1;
        end, clear i,

        histi = count;
        hConc(hStep,:) = histi'/dz*L/numel(pp);

    end, clear hStep pp histi,
    meanConc = mean(hConc, 'omitnan');
    stdConc = std(hConc, 'omitnan');
end