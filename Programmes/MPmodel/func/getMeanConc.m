function [meanConc, stdConc, hConc] = getMeanConc(zPart, nCat, dz)
% zPart = zFinal(end-10:end);
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

        histi = count;
        hConc(hStep,:) = histi'/dz;

    end, clear hStep pp histi,
    meanConc = mean(hConc, 'omitnan');
    stdConc = std(hConc, 'omitnan');
end