function ustar = ustarBeaufort(BeaufortForce, H)
    kappa = 0.40;   
    load(['/media/ian/Transcend/MPsDistrib/Data/Turbulence/DataTurbForce' num2str(BeaufortForce) '.mat'])
    Kmax = max(K);
    ustar = Kmax/(kappa*H)*4/3;
end