function saveTurb(Beaufort, H, N, Lat)
w_speed = 3*Beaufort.^(3/2);

z = 0:H/N:H;
z_ = (z(1:end-1)+z(2:end))/2;

date = datetime(2020,03,18);
[~,Row_day,~,z__day] = KsSalTemp2020(w_speed, date);
rhof = interp1(-z__day,Row_day,z,'pchip'); % density of sea water [kg.m⁻³]
clear KZ_day Row_day z_day z__day,

[K,dK,u_result,v_result]=Turbulence(H,N,rhof,Lat,w_speed)
save(['../Data/Turbulence/test_DataTurbForce' num2str(Beaufort) '.mat'])
end
