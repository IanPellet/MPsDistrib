function [MeanContent] = GetMeanContent(zHistory, iStart, iEnd, N, dz, AggContent)
    MeanContent = NaN(length(zHistory), N);
    for hStep = 1:length(zHistory)
        aCont = AggContent{hStep};
        aCont = aCont(iStart:iEnd);
        pp = zHistory{hStep};
        [ppSort, iSort] = sort(pp(iStart:iEnd));

        topBound = dz;
        j = 1;
        count = cell(N,1);

        for i=1:length(ppSort)
            while ppSort(i) > topBound
                topBound = topBound+dz;
                j = j+1;
            end
            count{j} = [count{j} aCont(iSort(i))];
        end, clear i,
        for ic = 1:N
            MeanContent(hStep,ic) = mean(count{ic});
        end
% if hStep >= 1000
% disp('pouet')
% end
    end
end