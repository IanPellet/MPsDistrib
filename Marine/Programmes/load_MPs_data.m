function [T] = load_MPs_data(data_file)
%LOAD_MPS_DATA Loads data file for the MPs
%   Loads a .txt file and sets the variable type of the columns.
    opts = detectImportOptions(data_file);
    opts = setvartype(opts,{'station'},'string');
    opts = setvartype(opts,{'replica', 'type', 'n'},'uint32');
    opts = setvartype(opts,{'date', 'time'},'datetime');
    disp([opts.VariableNames' opts.VariableTypes']);
    T = readtable(data_file, opts);
end

