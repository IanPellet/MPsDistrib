classdef MP < Particle
    %MP Micro-plastic class, represents one particle
    %
    properties (Constant)
        Type = 'MP' % Type of particle, constant = 'MP'
    end
    properties
        Locked = 0
    end
    
    methods
        function obj = lock(obj, agg)
            % Lock MP on OrgaAggr
            obj.Locked = agg.Index;
        end
    end
end

