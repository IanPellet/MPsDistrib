function [densite] = CalculDensiteMP(T,S)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %densite = 1000-0.12*T+0.35*S;
    Ro0 = 999.842594 ...
            + 6.793952E-2*T - 9.095290E-3 * T.^2 + 1.001685E-4 * T.^3 ...
            -1.120083E-6 * T.^4 + 6.536332E-9 * T.^5;
    A = 8.24493E-1 - 4.0899E-3 * T+7.6438E-5 * T.^2 - 8.2467E-7 * T.^3 ...
       + 5.3875E-9 .* T .^4;
    B = -5.72466E-3 + 1.0227E-4 * T - 1.6546E-6 * T.^2;
    C = 4.8314E-4;
    
    densite = Ro0 + A.*S + B.*S.^1.5 + C.*S.^2;    
    
end

