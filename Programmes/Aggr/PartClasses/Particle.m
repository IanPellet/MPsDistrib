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
    end
    properties (Dependent)
        Ws % Settling velocity (m.s⁻¹)
    end
    properties (Abstract, Constant)
        Type
    end

    methods
        %% Constructor
        function obj = Particle(size, rhop, rhow)
            %PARTICLE Constructor
            %   If no argument is passed (default constructor), everything
            %   is set to 0.
            %
            % PARAMETERS
            % size : double, particle size (m)
            % rhop : double, particle density (kg.m⁻³)
            % rhow : double, water density ((kg.m⁻³)
            
            if nargin == 0 % default constructor
                obj.Size = 0;
                obj.RhoP = 0;
                obj.RhoW = 0;
            else
                obj.Size = size;
                obj.RhoP = rhop;
                obj.RhoW = rhow;
            end
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


