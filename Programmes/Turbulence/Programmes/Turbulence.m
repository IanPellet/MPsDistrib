function [K,dK,u_result,v_result]=Turbulence(L,N,rhoW,Lat0,W_vent)

% Column discretisation
dz = L/N; 
z=0:dz:L; 
z_=(z(1:end-1)+z(2:end))/2; 

%Constantes
kappa=0.4; 
Kmol=10^(-6);
rho_air=1.292 ;

%Force de Coriolis
omega_T=7.2921*10^(-5);
f=2*omega_T*sin(Lat0);

%Coefficient de Trainée (surface)
z0=0.0002;   
Cd=(kappa/log(10/z0))^2;   

%Vitesse du vent en m/s
W=W_vent/3.6;


%%Systeme d'équation
syms u [1,N+1]
syms v [1, N+1]
syms K  [1,N]

%Initialisation et conditions aux limites
init_K=K(1)==Kmol+((kappa*z_(1)*(1-z_(1)/L))^2)*sqrt(((u(2)-u(1))/dz)^2+((v(2)-v(1))/dz)^2);
condition_surface1=K(1)*(u(2)-u(1))/dz==rho_air*Cd*W^2/rhoW(1);
condition_surface2=K(1)*(v(2)-v(1))/dz== 0;
condition_fond1=u(N+1)==0;
condition_fond2=v(N+1)==0;
equation=[condition_surface1,condition_surface2,init_K,condition_fond1,condition_fond2];
for i=2:N
eqn1=f*u(i)==(K(i)*(v(i+1)-v(i))/dz^2)-(K(i-1)*(v(i)-v(i-1))/dz^2);
eqn2=-f*v(i)==(K(i)*(u(i+1)-u(i))/dz^2)-(K(i-1)*(u(i)-u(i-1))/dz^2);
eqn3=K(i)==Kmol+((kappa*z_(i)*(1-z_(i)/L))^2)*sqrt(((u(i+1)-u(i))/dz)^2+((v(i+1)-v(i))/dz)^2);
equation=[equation,eqn1,eqn2,eqn3];
end 

%Resolution
Resolution=vpasolve(equation,[K,v,u]);
Result=structfun(@double,Resolution);

u_result=Result(N+1:2*N+1); 
v_result=Result(2*N+2:3*N+2);
%Calcul K et dK
K= Result(1:N);%Interpolation de K sur les extrémités
K2= interp1(z_,K,z,'makima');
dK=diff(K2)/dz;
%dK=diff(K);
%plot(K2,-z)
% hold on 
% xlabel('Viscosité turbulente (m^2/s')
% ylabel('profondeur(m)')
%plot(K,-z)
end