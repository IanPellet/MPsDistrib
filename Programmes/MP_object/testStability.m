function [dCmax, dC] = testStability(zPart, z, z_, dz, L)

    if length(zPart) < 4
        error("Not enought history time to test stability")
    end

    partition = ceil(length(zPart)/2);

    [meanConcTest, ~] = getMeanConc(zPart(1:partition), length(z_), dz);
    [meanConcFinal, ~] = getMeanConc(zPart(partition+1:end), length(z_), dz);

    dC = abs(meanConcTest-meanConcFinal)/mean(meanConcFinal);
    dCmax = max(dC);

    
end