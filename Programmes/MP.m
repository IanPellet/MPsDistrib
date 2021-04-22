classdef MP
    %MP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties     
        z_ % initial input, modified externaly
%         Ks_ % external profile
%         dKs_ % external profile
        row_
    end
    
    properties (SetAccess = private)
        type_ % input
        size_ % input
    end
    
    properties (Dependent)
        rho_ % dependent of type
        U_ % dep of type
        uz_ % dep of U and z
%         Kz_ % dependent of z and Ks
%         dKz_ % dependent of z and dKs
    end
    
    methods
        function obj = MP(type,size, z_init, row)
            
            if nargin == 0 % default constructor
                obj.type_ = 0;
                obj.size_ = 0;
                obj.z_ = 0;
                obj.row_ = 0;
            else
                obj.type_ = type;
                obj.size_ = size;
                obj.z_ = z_init;
                obj.row_ = row;
            end
            %MP Construct an instance of this class
            %   Size, type and initial position are given 
            %obj.Ks_ = Ks;
            %obj.dKs_ = dKs;
            
%             get.rho_(obj);
%             get.U_(obj);
%             
%             get.uz_(obj);
            %eval_Kz(obj);
            %eval_dKz(obj);
        end
        
        
        function rho = get.rho_(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
%             load('rho_per_type', 'RHO_TYPE'); % loads a variable RHO_TYPE 
            RHO_TYPE = [ 900 1050 ;
                        1000 1300 ;
                        1300 1400 ];
            % RHO_TYPE = [rho_min1 rho_max1 ;
            %             rho_min2 rho_max2 ;
            %             rho_min3 rho_max3 ]
            m = RHO_TYPE(obj.type_,1); % min Rho
            M = RHO_TYPE(obj.type_,2); % max Rho
            pd = makedist('Normal','mu',(m+M)/2 ,'sigma', (M-m)/4 ); % Normal law, centered on (m+M)/2, 95.45% of values between m and M
            rho = random(pd);
        end
        
        function U = get.U_(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Nom=[... 
                ;{'Nielsen'}...%Nielsen (1992)
                ;{'Soulsby'}...%Soulsby (1997)
                ;{'Ahrens'}...%Ahrens (2000)
               ];
            indNom = 3;
            S = obj.rho_./obj.row_; 
            Ws = eval(['Vitesse' cell2mat(Nom(indNom)) '(obj.size_,S);']);
            U = Ws; 
            test = obj.rho_<obj.row_;
            U(test) = -Ws(test);
        end
        
        function uz = get.uz_(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            uz = obj.U_(index);
        end
        
        
%         function obj = eval_Kz(obj)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             obj.Kz_ = obj.Ks_(obj.z_);
%         end
%         
%         function obj = eval_dKz(obj)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             obj.dKz_ = obj.dKs_(obj.z_);
%         end
    end
end

