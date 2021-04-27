function [dC] = DeltaConcentration(pos1,pos2,z,dz,dt)
%DELTACONCENTRATION Diff of concentration between two steps
% 

h1 = histogram(pos1, "BinEdges", z, 'Visible', 'off').Values;
h2 = histogram(pos2, "BinEdges", z, 'Visible', 'off').Values;
C1 = h1/dz;
C2 = h2/dz;

nPart = length(pos1);

dC = max(abs(C2 - C1)/dt)/nPart;
        
end

