function [znew] = eq3_NaiveRandomWalk(zi,Kzi,dt)
%EQ3_NAIVERANDOMWALK Diffusivity model
%   

    pd = makedist('Uniform', 'lower', -1, 'upper', 1);
    r = 1/3;
    R = random(pd,size(zi));
    znew = zi + R.*sqrt(2/r*Kzi*dt);
    
end

