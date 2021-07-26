function [new_x] = Step_Lagrangien_noDist(x, u, Kz, dKz, dt, L, pd)
%STEP_LAGRANGIEN One step foward in time for the Lagrangian transport model
%   The new value of x after dt is computed and returned by this function.
%disp(u)
%disp(Nu)
%disp(x)
% 
% pd = makedist('Uniform', 'lower', -1, 'upper', 1);
% d = 1/3;
% pd = makedist('Normal', 'mu', 0, 'sigma', d);
% d = 1;
d = pd.sigma;
R = random(pd,size(x));
new_x = x + dt*u + dKz*dt + R.*(2/d*Kz*dt).^(1/2);

% First-type boundary conditions
% new_x = max(0,new_x);
% new_x = min(L, new_x);

% Reflective boundary conditions
new_x(new_x < 0) = -new_x(new_x < 0);
new_x(new_x > L) = L-new_x(new_x > L) + L;

end



