filename = 'IO\store_INCAP_input_baseline.txt';
delimiter = '\t';
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
rawNumericColumns = raw(:, [1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]);
rawCellColumns = raw(:, 4);
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
storeINCAPinputbaseline = raw;
clearvars fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

%% corrections

% Correct Air Temp

% Correct Inflow temp

% Correct TP
% 
% %% write out
new_filename = 'IO\store_INCAP_input.txt';
copyfile(filename,new_filename)

fileID = fopen(new_filename,'w');
formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \r\n';
[nrows,ncols] = size(storeINCAPinputbaseline);
fprintf(fileID,formatSpec,storeINCAPinputbaseline{1,:}); % writing header
fprintf(fileID,formatSpec,storeINCAPinputbaseline{2,:}); % writing header

% because the global radiation column is full of 'NaN' that are stored as
% text ... 
temp=storeINCAPinputbaseline(3:end,4);
numind = cellfun(@isnumeric, temp);
temp(~numind) = {NaN};
storeINCAPinputbaseline(3:end,4)=temp;
input_write_out = (storeINCAPinputbaseline(3:end,:));
 
dlmcell(new_filename,input_write_out,'-a')
% 
