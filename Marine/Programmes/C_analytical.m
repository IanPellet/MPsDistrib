function [C_calc] = C_analytical(Ws, Ks, x_, CI)
%C_ANALYTICAL Returns the solution for C(z) at equilibrium
%   

C0 = sum(CI)/sum(exp(-Ws*x_/Ks));
C_calc = C0*exp(-Ws*x_/Ks);

end

