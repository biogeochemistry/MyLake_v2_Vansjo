%% Response surfaces input generator
clear all
tic

m_start=[1983, 1, 4]; % 1973 1 1 is the earliest possible given the provided input files8
m_stop=[2014, 5, 31]; % 2012, 12, 31 is the latest possible given the provided input files

num_days_2nd_half =(datenum(m_stop)-datenum(m_start))/2 ; % for the 2nd half of the run where we reduce TP inputs
days_frac_2nd_half = 100/(num_days_2nd_half);

%% config
total_runs = 0;
T_Air_loop = 0;
TP_loop = 0;
TP_red_loop = 0;
grid_size = 5.2; % give in percent (5.2 works with 8000 runs)
scaling_TP_red_vec = [1:num_days_2nd_half];
scaling_TP_red_mat = zeros((datenum(m_stop)-datenum(m_start)+1),8000);

%% air temperature T_Air
min_T_Air = -3;
max_T_Air = 10;

step_T_Air = (max_T_Air - min_T_Air) / (100/grid_size);

%% TP load
% this one is a log scale
min_TP = 0; % 0.01
max_TP = 3.3; % 
step_TP = (max_TP - min_TP) / (100/grid_size); % thus in log units

%% TP reduction
% this one is a log scale
max_TP_red = 10^max_TP/num_days_2nd_half/100; % this is the reduction per day needed max
step_TP_red = 2 * max_TP_red /(100/grid_size); % we will either increase or decrease

%% Opening  baseline input file
filename = 'IO\store_INCAP_input_baseline.txt';
delimiter = '\t';
startRow = 3;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
storeINCAPinputbaseline_mat= [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Opening header array
filename = 'IO\store_INCAP_input_baseline.txt';
delimiter = '\t';
endRow = 2;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow, 'Delimiter', delimiter, 'ReturnOnError', false);
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
storeINCAPinputbaseline_header = raw;
clearvars  delimiter endRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

for scaling_T_Air = min_T_Air:step_T_Air:max_T_Air
    T_Air_loop = T_Air_loop + 1;
    
    %correction
    storeINCAPinputbaseline_mat(:,6) = storeINCAPinputbaseline_mat(:,6) + scaling_T_Air;
    
    for new_TP_conc = min_TP:step_TP:max_TP
        TP_loop = TP_loop + 1;
        
        storeINCAPinputbaseline_mat(:,15) = 10^new_TP_conc;
        
       for scaling_TP_red = -max_TP_red:step_TP_red:max_TP_red % loop for the reduction or increase of TP in the 2nd half of the sim
            
            TP_red_loop = TP_red_loop + 1;
            
            % computing scaling vector for TP time-series in the 2nd half of the run
            scaling_TP_red_vec(1) = 10^new_TP_conc; % setting the initial TP from the outer loop
            for i = 2:length(scaling_TP_red_vec) % creating the rolling reduction vector ... sure it can be vectorized
                scaling_TP_red_vec(i) = scaling_TP_red_vec(i-1)-(scaling_TP_red_vec(1)*scaling_TP_red/100) ; % from the initial TP, creating the temporal trends
            end
            
            storeINCAPinputbaseline_mat(end-length(scaling_TP_red_vec)+1:end,15) = scaling_TP_red_vec;
            
            total_runs = total_runs + 1% run counter
            scaling_TP_red_mat(:,total_runs)=storeINCAPinputbaseline_mat(:,15); % stores the TP reduction for post-processing
          
          
            %% write out input file
%% uncomment for writing           
            new_filename = strcat('IO\RS\store_RS_input_',num2str(T_Air_loop),'_',num2str(TP_loop),'_',num2str(TP_red_loop),'.txt');
            list_input_file{total_runs,1} = new_filename;
            fileID = fopen(new_filename,'w');
            formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \r\n';
            fprintf(fileID,formatSpec,storeINCAPinputbaseline_header{1,:}); % writing header
            fprintf(fileID,formatSpec,storeINCAPinputbaseline_header{2,:}); % writing header
            dlmwrite(new_filename,storeINCAPinputbaseline_mat,'-append','delimiter','\t')
            fclose(fileID);


            
        end
        TP_red_loop = 0;
    end
    TP_loop = 0;
end



%% Running the model through all the input files
total_runs
clearvars -except list_input_file K_values* m_* scaling_TP_red_mat min_T_Air max_T_Air min_TP max_TP
toc