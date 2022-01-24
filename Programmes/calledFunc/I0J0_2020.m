function [IJ_tab] = I0J0_2020(LonLat_tab)
%% I0J0_2020 Get indices corresponding to lon-lat table (2020 data)
    ncFile = netcdf('../Data/rhoma2020/2020RHOMA_WRF6h_003.nc');
    Lon = ncFile{'longitude'}(:,:);
    Lat = ncFile{'latitude'}(:,:);
    
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