function [C_calc] = C_analytical(Ws, Ks, x_, N_init, L)
%C_ANALYTICAL Returns the solution for C(z) at equilibrium
%   

%C0 = sum(CI)/sum(exp(-Ws*x_/Ks));
C0 = Ws/Ks*N_init/(1-exp(-Ws*L/Ks));
%C0 = 0.194294734664083;
%disp(C0)
C_calc = C0*exp(-Ws*x_/Ks);

end

