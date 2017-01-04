function time_vec = fn_import_time_vec(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   VARNAME1 = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   VARNAME1 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   time_vec = fn_time_vec('volumeflowout.csv',1, 8214);
%

%% Initialize variables.
delimiter = '';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
time_vec = dataArray{:, 1};

