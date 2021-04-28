function [znew] = eq6_CorrectedRandomWalk(zi,Kzi,dKzi,K,dt,dz)
%EQ6_CORRECTEDRANDOMWALK Diffusivity model
%   

    pd = makedist('Uniform', 'lower', -1, 'upper', 1);
    r = 1/3;
    R = random(pd,size(zi));
%     ziOffset = zi + dKzi*dt/2;
%     index = max(1, fix(ziOffset/dz));
%     KziOffset = K(index);
%     znew = zi + dKzi*dt + R.*sqrt(2/r*KziOffset*dt);
    znew = zi + dKzi*dt + R.*sqrt(2/r*Kzi*dt);
    
end

