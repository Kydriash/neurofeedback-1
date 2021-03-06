function [data, col_names, row_names, first_col_name] = read_table(fpath_in)
% Reads table: 1-st column - lastnames, 1-st row - names of variables
% USAGE: [data, col_names, row_names, first_col_name] = read_table(fpath_in)
% Input:
% fpath_in - path to the tab-delimited text file containing table
% Output:
% data - numeric matrix with dable inner data
% col_names - names of columns (from 1-st row)
% row_names - names of rows (from 1-st column)
% first_col_name - value in the upper-left field


% Open file with input data
fid = fopen(fpath_in, 'r');

% Read header line
Q = textscan(fid, '%s', 1, 'delimiter', '\n', 'bufsize', 16000);

% Parse header line: get variable names
col_names = regexpi(Q{1}{1}, '\t', 'split');
first_col_name = col_names{1};
col_names = col_names(2:end);
ncols = length(col_names);

% Template for parsing line of data
templ = ['%s\t' repmat('%f\t', 1, ncols-1) '%f\n'];

row_names ={};
data = [];

while ~feof(fid)
    
    % Read line of data
    %Q = textscan(fid, templ, 'TreatAsEmpty', '\t');
    Q = textscan(fid, templ, 'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
    
    
    % Get row name
    row_names = [row_names Q{1}];
    
    % Get numeric data
    X = Q(2:end);
    X = [X{:}];
    data = [data; X];    
    
end

% Close input data file
fclose(fid);


end

