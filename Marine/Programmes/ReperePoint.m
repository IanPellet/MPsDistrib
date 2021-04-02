function [I0,J0]=ReperePoint(Lon,Lat,Lon0,Lat0)

%distance = Distance(Lon,Lat,Lon0,Lat0);
%[dd,II]=min(distance);
%[dd,J0]=min(min(distance));
%I0=II(J0);
%d0=distance(I0,J0);

d_lon = abs(Lon-Lon0);
d_lat = abs(Lat-Lat0);
[~, J0] = min(d_lon(1,:));
[~, I0] = min(d_lat(:,1));

end
