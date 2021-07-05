function [new_x, new_y, new_z] = Step_Aggr(x, y, z, u, Kz, dKz, dt, L)
%STEP_LAGRANGIEN One step foward in time for the Lagrangian transport model
%   The new value of x after dt is computed and returned by this function.
%disp(u)
%disp(Nu)
%disp(x)
% 
pd = makedist('Uniform', 'lower', -1, 'upper', 1);
d = 1/3;
% pd = makedist('Normal', 'mu', 0, 'sigma', 1);
% d = 1;

R = random(pd,size(x));
new_x = x + (dKz*dt + R.*(2/d*Kz*dt).^(1/2))/L;
new_x(new_x < 0) = -new_x(new_x < 0);
new_x(new_x > 1) = 1-new_x(new_x > 1) + 1;

R = random(pd,size(y));
new_y = y + (dKz*dt + R.*(2/d*Kz*dt).^(1/2))/L;
new_y(new_y < 0) = -new_y(new_y < 0);
new_y(new_y > 1) = 1-new_y(new_y > 1) + 1;


R = random(pd,size(z));
new_z = z + dt*u + dKz*dt + R.*(2/d*Kz*dt).^(1/2);

% First-type boundary conditions
% new_x = max(0,new_x);
% new_x = min(L, new_x);

% Reflective boundary conditions
new_z(new_z < 0) = -new_z(new_z < 0);
new_z(new_z > L) = L-new_z(new_z > L) + L;

end
