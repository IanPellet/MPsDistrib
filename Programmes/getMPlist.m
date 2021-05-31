function [mp_list] = getMPlist(nPart, sizePart, rhop, rhow, rFrag)
    % allocate memory to store particles
    mp_list(nPart) = MP; % array of MP objects
    % Fill the array
    for i = 1:(nPart)
        mp_list(i) = MP(sizePart(i), rhop, rhow, rFrag);
    end, clear i,
end