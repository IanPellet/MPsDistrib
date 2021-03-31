function [new_x] = Step_Lagrangien(u, Nu, x)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global dt

new_x = x + dt*u + Nu*dt.^2;
new_x = max(0,new_x);

end

