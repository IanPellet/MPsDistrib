classdef MP
    %MP Micro-plastic class, represents one particle
    %  
    properties
        rhow_
        size_ % (m)
        sticky_ = 0.5
        aggr_ = false
        dp_ = 0
    end
    properties (SetAccess = private)
        type_ % int {0,1,2,3} NON IMPLÉMENTÉ
        
        rho_ % (kg.m⁻³)
        U_ % velocity in the water column, dep of type (m.s⁻¹)
        fragRate_ % double in (0,1], probability for a particle to be fragmented each second
    end

    methods
        %% Constructor
        function obj = MP(size, rhop, rhow, fragRate)
            %MP Constructor
            % If no argument is passed (default constructor), everything
            % is set to 0.
            %
            % PARAMETERS
            % size : double, particle size (m)
            % rhop : double, particle density (kg.m⁻³)
            % rhow : double, water density ((kg.m⁻³)
            % fragRate : double in (0,1], probability for a particle to be fragmented
            % each second 
            %
            
            obj.type_ = 0; % TYPE NOT IMPLEMENTED YET
            
            if nargin <= 1 % default constructor
                obj.size_ = 0;
                obj.rho_ = 0;
                obj.rhow_ = 0;
                obj.fragRate_ = 0;
            else
                obj.size_ = size;
                obj.rho_ = rhop;
                obj.rhow_ = rhow;
                obj.fragRate_ = fragRate;
            end           
        end
        
        %% Setter
        function value = get.U_(obj)
            %SETU Compute fall velocity of particle in the water column
            %   
            Nom=[... 
                ;{'Nielsen'}...%Nielsen (1992)
                ;{'Soulsby'}...%Soulsby (1997)
                ;{'Ahrens'}...%Ahrens (2000)
               ];
            indNom = 3;
            S = obj.rho_./obj.rhow_; 
            Ws = eval(['Vitesse' cell2mat(Nom(indNom)) '(obj.size_,S);']);
            value = Ws'; 
            test = obj.rho_<obj.rhow_;
            value(test) = -Ws(test);
            
            if obj.aggr_
                value = VitesseLiLogan(obj.size_, obj.dp_, 2.5, value);
            end
        end

    end
end

