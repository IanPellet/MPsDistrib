classdef Aggr < Particle
    %AGGR Organic aggregate particle
    %   Particle subclass representing organic mater aggregates that can
    %   trap MP particles
    properties (Constant)
        Type = 'Aggr' % Type of particle, constant = 'AggrOrga'
    end
    properties
        Adh % particle adherence
    end
    properties (SetAccess = private)
        Content % MP contained by the aggregate
    end
    methods
        %% Constructor
        function obj = Aggr(size, rhop, rhow, index, adherence)
            %AGGR Constructor
            %   If no argument is passed (default constructor), everything
            %   is set to 0.
            %
            % PARAMETERS
            % size : double, particle size (m)
            % rhop : double, particle density (kg.m⁻³)
            % rhow : double array, water density (kg.m⁻³)
            % index : int, particle index (to identify it in a simulation)
            % adherence : double, [0,1] probability to stick when touching
            % other particle
            
            if nargin == 0 % default constructor
                size = 0;
                rhop = 0;
                rhow = 0;
                index = 0;
                adherence = 0;
            end
            
            % Call Particle constructor
            obj@Particle(size, rhop, rhow, index);
            
            obj.Adh = adherence;
            obj.Content = [];
        end
        
        function [obj, mp] = aggrMP(obj, mp)
            % AGGRMP aggregation of one MP particle into the aggregate
            %
            % PARAMETERS
            % obj : Aggr object
            % mp : MP object
            
            mp = mp.lock(obj); % set MP as locked on the Aggr
            obj.Content = [obj.Content mp]; % add MP to the content of the Aggr
        end
    end
end


