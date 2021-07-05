function [K, dK] = turbulence(z, dz)

    K = k(z);
    dK = (k(z+dz/2)-k(z-dz/2))/dz;
    
    function k = k(z)
        K0 = 0.08;
        KWM = 1e-5;
        zWM = 20; % depth of the well mixed layer (m)
        k(z < zWM) = K0 + z(z < zWM).*(KWM-K0)./zWM;
        k(z >= zWM) = KWM;    
    end
end