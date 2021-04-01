function [new_x] = Step_Lagrangien(x, u, Nu)
%STEP_LAGRANGIEN One step foward in time for the Lagrangian transport model
%   The new value of x after dt is computed and returned by this function.
global dt
%disp(u)
%disp(Nu)
%disp(x)

new_x = x + dt*u ;
new_x = max(0,new_x);

end

