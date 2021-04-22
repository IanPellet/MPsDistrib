classdef WC
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        
        
    end
    
    methods
        function obj = WC(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            %% Load hydrodinamic model
            ModeleHydro='2012RHOMA_arome_003.nc';
            SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
            load(SauvegardeModeleHydro)
            
            %% Water column parameters
            L = 50;
            N=50;  dz= L/N;  z=0:dz:L; % z : boundaries of the meshes
            z_=(z(1:end-1)+z(2:end))/2; % middle of each mesh
            z0=L*Sigma;

            %% Turbidity initialisation
            DensiteFevrierRhoma 
            [K,dK] = wcp_interpolation(z0,KZ_Fev10,-z_); % Diffusivity
            % K_val = 0.01;
            % K = ones(size(z))*K_val;
            % dK = zeros(size(z));
            % K = K*cKs;

        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

