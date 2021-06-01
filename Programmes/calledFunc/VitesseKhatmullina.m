function [ w ] = VitesseKhatmullina( d,~,~, L ,g_red)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

w = pi/2 * g_red * d*L/(55.235*L+12.691);

end