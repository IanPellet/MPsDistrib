% Creat video of moving particles

BeaufortForce = 2;

% v=3.B^{3/2} si v  est en km/h, Wikipedia
wind_speed = 3*BeaufortForce.^(3/2);
H = 70;
tf = 60*60*24*3; 
dt = 60;
P = 15000; % Number of particles
dt_test = 60*60;
partZinit = linspace(0, H, P);

%% create list of particles
Prho = 1000;
rhoB = 1048;
Psize = 750e-6;
fb = 0;
% allocate memory to store particles
part_list(P) = Particle; % array of Particle objects
% Fill the array
L0 = min(Psize);
for iP = 1:(P)
    part_list(iP) = Particle(Psize, rhoB, Prho, fb, iP, L0); % size and densities only used if ws computed
end, clear iP,

load(['../Data/Turbulence/DataTurbForce' num2str(BeaufortForce) '.mat'])

date = datetime(2020,03,18);
[~,Row_day,~,z__day] = KsSalTemp2020(wind_speed, date);
rhof = interp1(-z__day,Row_day,z,'pchip'); % density of sea water [kg.m⁻³]
clear Row_day z__day date,

posPart = Part_Simulator(part_list, partZinit, K, dK, H,...
    dz, rhof, tf, dt_test, dt, 1e-10, true, false, false, false, 0);


%% 
Conc = NaN(length(posPart),N);
hStep = 1;

    pp = posPart{hStep};
    ppSort = sort(pp(1:length(pp)));

    topBound = dz;
    j = 1;
    count = zeros(N,1);

    for i=1:length(ppSort)
        while ppSort(i) > topBound
            topBound = topBound+dz;
            j = j+1;
        end
        count(j) = count(j)+1;
    end, clear i,
    Conc(hStep,:) = count/dz;

plot(Conc,-z_)
ylim([-H 0])
xlabel("Concentration (particles.m⁻¹)")
ylabel("Depth (m)")
