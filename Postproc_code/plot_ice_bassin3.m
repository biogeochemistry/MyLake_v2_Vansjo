%% RMC June 2013 rmc@niva.no
init_datum = datenum(1973,1,1);

%% import MyLake output
% filename = 'basin1freezingdates.csv';
% delimiter = '';
% formatSpec = '%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% Basin1_on = dataArray{:, 1}+init_datum;

% filename = 'basin2freezingdates.csv';
% delimiter = '';
% formatSpec = '%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% Basin2_on = dataArray{:, 1}+init_datum;

filename = 'basin3freezingdates.csv';
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
Basin3_on = dataArray{:, 1}+init_datum;

% filename = 'basin1meltingdates.csv';
% delimiter = '';
% formatSpec = '%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% Basin1_off = dataArray{:, 1}+init_datum;

% filename = 'basin2meltingdates.csv';
% delimiter = '';
% formatSpec = '%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% Basin2_off = dataArray{:, 1}+init_datum;
% 
filename = 'basin3meltingdates.csv';
delimiter = '';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
Basin3_off = dataArray{:, 1}+init_datum;

clearvars filename delimiter formatSpec fileID dataArray ans;

%% Create sigma-plot cols yyyy on off duration

Basin3 = zeros(50, 4);
maxo = length(Basin3_on);
maxf = length(Basin3_off);
Basin3(1:maxo,1) = str2num(datestr(Basin3_on(:,1),'YYYY')); % extracts year from datenum
temp_year= datenum(Basin3(1:maxo,1),1,1) ; % datenum the 1st day of the year
Basin3(1:maxo,2) = Basin3_on(:,1) - temp_year(1:maxo,1); % creates ice-on dates
Basin3(1:maxf,3) = Basin3_off(:,1) - temp_year(1:maxf,1); % creates ice-off dates
Basin3(:,4) = Basin3(:,3) - Basin3(:,2); % creates duration
Basin3(1,:) = []; % delete duplicate 1st yr

% checks double years and fall ice-off events
for a = 1:maxo-1; 
    if Basin3(a,1)==Basin3(a+1,1) % if the next year is same as current
        Basin3(a,3)=Basin3(a+1,3); % replace current yr ice-off by next
        Basin3(a+1,:)=[]; % delete second year 
        Basin3(maxo,:)=0; % make sure matrix keeps size
    end    
end

    


%% Other basins


%% 
%% Calculate abs mean deviation
% Please prepare observation matrix ... yyyy, on, off, (duration). Days are
% elapsed since Jan 1st of the year Ice went on. 
load ice_obs.mat;
% performance for basin1
compare = intersect(Basin3(:,1),ice_obs(:,1)); %create vector with yr common to model and obs
j = 0;
for i = 1:length(compare) % compare contains those yr that are common to both sets
    index_m = Basin3(:,1)==compare(i,1); % finds the ith comparable yr in model
    index_s = ice_obs(:,1)==compare(i,1); % finds the ith comparable yr in obs
    if ice_obs(index_s,2) ~= 0 % will not calculate dev when the is no obs
        j = j+1;
        dev3_on(j,1) = ice_obs(index_s,1);
        dev3_on(j,2) = Basin3(index_m,2)-ice_obs(index_s,2); % stores dev ice_on
    end
end
dev3_on(j+1,2)= mean(dev3_on(:,2));
j = 0;
for i = 1:length(compare) % compare contains those yr that are common to both sets
    index_m = Basin3(:,1)==compare(i,1); % finds the ith comparable yr in model
    index_s = ice_obs(:,1)==compare(i,1); % finds the ith comparable yr in obs
    if ice_obs(index_s,3) ~= 0 && Basin3(index_m,3) ~= 0;  % will not calculate dev when the is no obs
        j = j+1;
        dev3_off(j,1) = ice_obs(index_s,1);
        dev3_off(j,2) = Basin3(index_m,3)-ice_obs(index_s,3); % stores dev ice_on
    end
end

dev3_off(j+1,2)= mean(dev3_off(:,2));
    
clear Basin3_on Basin3_off compare a i index_m index_s maxo maxf temp_year init_datum j
