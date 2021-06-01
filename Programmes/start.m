% add curent folder and its subdirectories 
currDir = genpath(pwd);
addpath(currDir) 

% add path to data
addpath('../Data/')

% save added paths
savepath pathdef.m