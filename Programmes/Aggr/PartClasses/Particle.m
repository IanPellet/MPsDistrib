classdef (Abstract) Particle
    %PARTICLE Water column particle class
    %   Abstract class representing any kind of particle present in the 
    %   water column.
    
    properties (SetAccess = private)    
        Size % particle size (m)
        RhoP % particle density (kg.m⁻³)
    end
    properties
        RhoW % water density (kg.m⁻³)
        Index % list index of the particle
    end
    properties (Dependent)
        Ws % Settling velocity (m.s⁻¹)
    end
    properties (Abstract, Constant)
        Type % MP or aggregate
    end

    methods
        %% Constructor
        function obj = Particle(size, rhop, rhow, index)
            %PARTICLE Constructor
            %   If no argument is passed (default constructor), everything
            %   is set to 0.
            %
            % PARAMETERS
            % size : double, particle size (m)
            % rhop : double, particle density (kg.m⁻³)
            % rhow : double, water density (kg.m⁻³)
            % index : int, particle index (to identify it in a simulation)
            
            if nargin == 0 % default constructor
                size = 0;
                rhop = 0;
                rhow = 0;
                index = 0;
            end
            obj.Size = size;
            obj.RhoP = rhop;
            obj.RhoW = rhow;
            obj.Index = index;
        end
        
        %% Setter
        function value = get.Ws(obj)
            %SETWS Compute settling velocity of the particle in the water column
            %   Ahrens2000 settling velocity formula

            S = obj.RhoP./obj.RhoW; 
            value = VitesseAhrens(obj.Size,S)';
            test = obj.RhoP < obj.RhoW;
            value(test) = -value(test);           
            
        end

    end
end


