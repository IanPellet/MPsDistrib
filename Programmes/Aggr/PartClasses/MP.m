classdef MP < Particle
    %MP Micro-plastic class, represents one particle
    %
    properties (Constant)
        Type = 'MP' % Type of particle, constant = 'MP'
    end
    properties
        % when a MP particle is free Locked = 0, when it aggregates with
        % other particles, it's set to the index of the Aggr it's in
        Locked = 0 
    end
    
    methods
        function obj = lock(obj, agg)
            %LOCK Lock MP on Aggr
            obj.Locked = agg.Index;
        end
    end
end

