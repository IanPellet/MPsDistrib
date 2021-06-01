function [T] = load_ConcVol_data(data_file)
%LOAD_MPS_DATA Loads data file for the MPs
%   Loads a .txt file and sets the variable type of the columns.
    opts = detectImportOptions(data_file);
    opts = setvartype(opts,{'date'},'datetime');
    opts = setvartype(opts,{'station'},'string');
    T = readtable(data_file, opts);    
end

