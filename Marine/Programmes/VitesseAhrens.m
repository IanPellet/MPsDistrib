function [ w ] = VitesseAhrens( d, S, D_)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global g nuw 

Delta=abs(S-1);

A=(Delta.*g.*d.^3)/nuw^2;

Cl=0.055*tanh((12*A.^-0.59).*exp(-0.0004*A));
Ct=1.06*tanh((0.016*A.^0.5).*exp(-120./A));

w= (Cl.*Delta.*g.*d.^2)./nuw+Ct.*(Delta.*g.*d).^0.5;

end