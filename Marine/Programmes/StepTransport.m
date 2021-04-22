function [C,outputArg2] = StepTransport(u,Nu,C,schema);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global dt dz
if strcmp(schema,'CentreDF') 
    UC=NaN(1,size(C,2)+2);
elseif strcmp(schema,'Ordre4')
    UC=NaN(1,size(C,2)+4);
elseif  strcmp(schema,'UpWind') | strcmp(schema,'CentreVF')
    C_=NaN(1,size(C,2)+1);
end
Adv=0*C;Diff=0*C;

if strcmp(schema,'Ordre4')
    UC(4:end-3)=(u(1:end-1)+u(2:end)).*C(2:end-1)/2;
    UC(3)=u(1)*C(1)/2;     UC(2)=0;    UC(1)=0;
    UC(end-2)=u(end)*C(end)/2;  UC(end-1)=0;  UC(end)=0;   
    UC(3)=0;     UC(2)=0;    UC(1)=0;
    UC(end-2)=0;  UC(end-1)=0;  UC(end)=0;   
%         Adv(2) = dt/(dz)*(C(1) - C(2))+C(2);
%         Adv(3:end-2) = (dt)/(12*dz)*(u(5:end).*C(5:end) - 8*u(4:end-1).*C(4:end-1)...
%                                + 8*u(2:end-3).*C(2:end-3)- u(1:end-4).*C(1:end-4))...
%                                + C(3:end-2);
%         Adv(end-1) = dt/(dz)*(C(end-2) - C(end-1))+C(end-1);
    Adv = (dt)/(12*dz)*(UC(5:end) - 8*UC(4:end-1)...
                           + 8*UC(2:end-3)- UC(1:end-4));
    UC(2:end-1)=u.*C(1:end-1);
        ii=find(u<0);  UC(ii+1)=u(ii).*C(ii+1);
        UC(1)=0;    UC(end)=0;
    Adv(1:2) = dt/(dz).*(UC(1:2) - UC(2:3));
    Adv = dt/(dz).*(UC(1:end-1) - UC(2:end));
elseif  strcmp(schema,'CentreDF')
%   Cf=(C(1:end-1)+C(2:end))/2
    UC(3:end-2)=(u(1:end-1)+u(2:end)).*C(2:end-1)/2;
    UC(2)=u(1)*C(1)/2;    UC(1)=0; 
    UC(end-1)=u(end)*C(end)/2;  UC(end)=0;   
    %Adv(2:end-1) = dt/(2*dz)*(u(1:end-2).*C(1:end-2)-u(3:end).*C(3:end))+C(2:end-1);
    Adv = dt/(dz).*(UC(1:end-2) - UC(3:end));
elseif  strcmp(schema,'CentreVF')
    C_(2:end-1)=(C(1:end-1)+C(2:end))/2;
    C_(1)=0;C_(end)=0;
    Adv = dt/(dz).*(u(1:end-1).*C_(1:end-1) - u(2:end).*C_(2:end));
elseif  strcmp(schema,'UpWind')
    C_(2:end-1)=C(1:end-1);
    ii=find(u(1:end-1)<0);  C_(ii)=C(ii);
    C_(1)=0;C_(end)=0;
    Adv = dt/(dz).*(u(1:end-1).*C_(1:end-1) - u(2:end).*C_(2:end));
    %Adv = dt/(dz).*(u*C_(1:end-1) - u*C_(2:end));
end


Diff(2:end-1)=dt/(dz.^2) * (Nu(2:end-1).*C(3:end)...
                        - (Nu(2:end-1)+Nu(1:end-2)).*C(2:end-1)...
                        +Nu(1:end-2).*C(1:end-2));
Diff(1)=dt/(dz.^2) * Nu(1).*(C(2)- C(1)); 
Diff(end)=-dt/(dz.^2) * Nu(end).*(C(end)- C(end-1)); 

C = C + Adv + Diff;
end

