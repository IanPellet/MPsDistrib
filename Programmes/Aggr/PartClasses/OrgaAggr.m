classdef OrgaAggr < Particle
    %AGGRORGA Organic aggregate particle
    %   Particle subclass representing organic mater aggregates that can
    %   trap MP particles
    properties (Constant)
        Type = 'OrgaAggr' % Type of particle, constant = 'AggrOrga'
    end
    properties
        Adh % particle adherence
    end
    properties (SetAccess = private)
        Content % MP contained by the aggregate
    end
    methods
        %% Constructor
        function obj = OrgaAggr(size, rhop, rhow, adherence)
            %AGGRORGA Constructor
            %   If no argument is passed (default constructor), everything
            %   is set to 0.
            %
            % PARAMETERS
            % adherence : 
            
            if nargin == 0 % default constructor
                size = 0;
                rhop = 0;
                rhow = 0;
                adherence = 0;
            end
            
            % Call Particle constructor
            obj@Particle(size, rhop, rhow);
            
            obj.Adh = adherence;
            obj.Content = [];
        end
        
        function obj = aggrMP(obj, mp)
            % AGGRMP aggregation of one MP particle into the aggregate
            obj.Content = [obj.Content mp];
        end
    end
end


