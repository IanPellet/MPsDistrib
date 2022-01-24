function [IJ_tab] = I0J0(LonLat_tab)
%% I0J0 Get indices corresponding to lon-lat table (2012 data)
    ModeleHydro='2012RHOMA_arome_003.nc';
    SauvegardeModeleHydro=['DonneeBase' ModeleHydro(1:end-3)];
    load(SauvegardeModeleHydro)
    
    IJ_tab = LonLat_tab;
    IJ_tab.Properties.VariableNames = [{'Station'},{'I0'},{'J0'}];
    for row = 1:size(LonLat_tab,1)
        Lon0 = LonLat_tab(row,2).Variables;
        Lat0 = LonLat_tab(row,3).Variables;
        [I0,J0]=ReperePoint(Lon,Lat,Lon0,Lat0);
        IJ_tab(row,2).Variables = I0;
        IJ_tab(row,3).Variables = J0;
    end
end