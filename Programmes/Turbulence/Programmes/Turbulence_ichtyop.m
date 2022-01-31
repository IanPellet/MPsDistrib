function[K,dK]=Turbulence_ichtyop(L,u,v,z_,z,dz)
%u et v les vitesses d'advection horizontale, z_ les coordonnées des milieux
%des mailles, z aux extrémités des mailles et L la profondeur totale

%Constantes
Kmol=10^(-6);%viscosité moléculaire
kappa=0.4; %constante de karman

K=zeros(1,length(z_));
for i=1:length(z_)
K(i)=Kmol+((kappa*z_(i)*(1-z_(i)/L))^2)*sqrt(((u(i+1)-u(i))/dz)^2+((v(i+1)-v(i))/dz)^2);
end 

%Interpolation
K=interp1(z_,K,z,'pchip');
dK =diff(K);
end 