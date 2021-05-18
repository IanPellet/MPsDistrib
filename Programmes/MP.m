classdef MP
    %MP Micro-plastic class, represents one particle
    %   Detailed explanation goes here
    
    properties     
        z_ % particle's position (m)

        dz_
    end
    
    properties (SetAccess = private)
        type_ % int {0,1,2,3} NON IMPLÉMENTÉ
        size_ % (m)
        rho_ % (kg.m⁻³)
        U_ % velocity in the water column, dep of type (m.s⁻¹)
    end
    
    properties (Dependent)
        
        uz_
        index_
    end
    
    methods
        function obj = MP(size, rhop, rhow)
            %MP Constructor
            % If no argument is passed (default constructor), everything
            % is set to 0.
            %
            % PARAMETERS
            % size : double, particle size (m)
            % rhop : double, particle density (kg.m⁻³)
            % z : double, particle's position (m)
            % rhow : double, water density ((kg.m⁻³)
            %
            
            obj.type_ = 0; % TYPE NOT IMPLEMENTED YET
            
            if nargin == 0 % default constructor
                obj.size_ = 0;
                obj.rho_ = 0;
%                 obj.z_ = 0;
%                 obj.dz_ = 0;
                obj.U_ = 0;
            else
                obj.size_ = size;
                obj.rho_ = rhop;
%                 obj.z_ = z;
%                 obj.dz_ = dz;
                obj.U_ = obj.setU(rhow);
            end
            
            
            
        end
        

        function value = setU(obj,rhow)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Nom=[... 
                ;{'Nielsen'}...%Nielsen (1992)
                ;{'Soulsby'}...%Soulsby (1997)
                ;{'Ahrens'}...%Ahrens (2000)
               ];
            indNom = 3;
            S = obj.rho_./rhow; 
            Ws = eval(['Vitesse' cell2mat(Nom(indNom)) '(obj.size_,S);']);
            value = Ws'; 
            test = obj.rho_<rhow;
            value(test) = -Ws(test);
        end
        
%         function uzVal = get.uz_(obj)
%             uzVal = obj.U_(obj.index_);
%         end
%         
%         function iVal = get.index_(obj)
%             iVal = max(1, cast(obj.z_/obj.dz_, 'uint32'));
%         end
%         
%         function new_z = Step_Lagrangien(obj, K, dK, dt, L)
%             %STEP_LAGRANGIEN One step foward in time for the Lagrangian transport model
%             %   The new value of x after dt is computed and returned by this function.
%             
%             Kz = K([obj.index_]);
%             dKz = dK([obj.index_]);
%             
%             pd = makedist('Uniform', 'lower', -1, 'upper', 1);
%             d = 1/3;
%             R = random(pd,size(obj));
% 
%             new_z = obj.z_ + dt*obj.uz_ + dKz*dt + R.*(2/d*Kz*dt).^(1/2);
%             % First-type boundary conditions
%             % new_x = max(0,new_x);
%             % new_x = min(L, new_x);
% 
%             % Reflective boundary conditions
%             new_z(new_z < 0) = -new_z(new_z < 0);
%             new_z(new_z > L) = L-new_z(new_z > L) + L;
% 
%             end
        
        %         function rho = get.rho_(obj)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
% %             load('rho_per_type', 'RHO_TYPE'); % loads a variable RHO_TYPE 
%             RHO_TYPE = [ 900 1050 ;
%                         1000 1300 ;
%                         1300 1400 ];
%             % RHO_TYPE = [rho_min1 rho_max1 ;
%             %             rho_min2 rho_max2 ;
%             %             rho_min3 rho_max3 ]
%             m = RHO_TYPE(obj.type_,1); % min Rho
%             M = RHO_TYPE(obj.type_,2); % max Rho
%             pd = makedist('Normal','mu',(m+M)/2 ,'sigma', (M-m)/4 ); % Normal law, centered on (m+M)/2, 95.45% of values between m and M
%             rho = random(pd);
%         end
        
    end
end

