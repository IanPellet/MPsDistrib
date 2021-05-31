function [dCmax, dC] = testStability(zPart, z, z_, dz, L)

    if length(zPart) < 4
        error("Not enought history time to test stability")
    end

    partition = ceil(length(zPart)/2);

    [meanConcTest, ~] = getMeanConc(zPart(1:partition), z, z_, dz, L);
    [meanConcFinal, ~] = getMeanConc(zPart(partition+1:end), z, z_, dz, L);

    dC = abs(meanConcTest-meanConcFinal)/mean(meanConcFinal);
    dCmax = max(dC);

    
end