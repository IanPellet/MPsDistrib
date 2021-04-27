function [Part] = UpdatePart(Part,K,dK,dz)
%UPDATEPART Updates K and dK values corresponding to the part's position
%   
z_part = Part(1,:);
index = max(1, fix(z_part/dz));
K_part = K(index);
dK_part = dK(index);
Part = [z_part ; K_part ; dK_part];
end

