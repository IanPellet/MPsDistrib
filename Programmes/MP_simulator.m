function [zFinal] = MP_simulator(mp, zInit, K, dK, L, dz, tf, dt_test)
%MP_SIMULATOR Run simulation with particles in mp
%
% mp : MP objects array, modeled particles
% zInit : double array, initial particle's position (m)
% K : double array, diffusivity profile (m².s⁻¹)
% dK : double array, diffusivity gradient (m.s⁻¹)
% L : double, water column depth (m)
% dz : double, water column discretisation step (m)
% tf : double, simulation time (s)
% dt_test : double, test interval (s)
%
% zFinal : double array, final particle's position (m)
%

    U = [mp.U_]; % 2D double array, mp fall velocities on the column (m.s⁻¹)
    uz = NaN(size(mp)); % double array, mp fall velocities (m.s⁻¹)
    zPart = zInit; % double array, particle's position (m)
    
    %% Time step initialisation
    dt = 10; % double, time step (s)
    ddK = diff(dK)./dz; % double array, diffusivity gradient's derivative (s⁻¹)
    dt = min(dt, abs(min(1./ddK)/10)); % check condition dt<<min(1/ddK) 
    dt = min(dt, dz/max(max(abs(U)))); % check condition dt < dz/max|u|
    
    %% Simulation
    t=0; 
    OnContinue=true;
    while OnContinue
        % Time update
        t=t+dt;
    
        %% Particules update
        index = max(1, cast(zPart/dz, 'uint32')); % int array, index of each particle's current mesh
        % Find current fall velocity of each particle
        for i=1:length(zPart)
            uz = U(index(i),i);
        end
        % Update particle's postion
        zPart = Step_Lagrangien(zPart, uz, K(index), dK(index), dt, L);
    
        % Test
        if (mod(t,dt_test)<=dt/2 || dt_test-mod(t,dt_test)<=dt/2 )
            % Test if final time is reached
            if (t>=tf)
               OnContinue = false;
            end
            % Display progress
            if mod(t,60)<dt/2
                disp([' Temps : ' num2str(t/tf*100) '%'])
            end
        end

    end
    zFinal = zPart; % double array, final particle's position (m)
    
end