function [ w ] = VitesseKhatmullina( d,~,~, L ,g_red)
%VITESSEKHATMULLINA Khatmullina fall velocity


w = pi/2 * g_red * d*L/(55.235*L+12.691);

end