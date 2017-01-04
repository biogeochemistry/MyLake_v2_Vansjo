function [store_INCAP_input, vanem_INCAP_input, INCA_QC] = fn_INCA_MyL(run_INCA, run_ID, clim_ID)

% This script launches INCA, write a MyLake input from INCA, then launches MyLake. In a nutchell, this prepares the "land inputs" for MyLake

%% Run INCA if INCA_run selected
Working_folder = strcat('C:\OneDrive\PROJECT\O-12111_MARS\4.4_Northern_basin\Vansjø\Model\INCA_to_MyLake\INCA_P\Scenario_Runs\',run_ID,'\');

if run_INCA == 1 
    cd (Working_folder); % ##### must replace this by the directory list
   !INCA.bat
else
    cd (Working_folder);
    disp 'Using existing output'
end

disp 'reading INCA output ...'

%%  Load Reach 
fileID = fopen('inca_out.dsd');% We'll use MATLAB's textscan to read blocks with Reach#{(Q)(SS)(TDP)(PP)(TP)(SRP)(Temp)}
%             Q                                          SS     TDP PP                        TP SRP T
formatSpec = '%f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %*f %f %f %*f %*f %*f %*f %*f %*f %f %f %f %*f %*f %*f %*f'; % %*f ignors the float column

% the first reach is special, it has 8 lines before start
Reach{1,:} = textscan(fileID,formatSpec,'HeaderLines',8,'Delimiter','\t'); % takes the first "block" --> reach 1

%% Bug: date lenghts not the same

k = 1;
while ~feof(fileID) % all the next one just need to skip 4 lines, so we look until there is no more reaches
    k = k+1;
    Reach{k,:} = textscan(fileID,formatSpec,'HeaderLines',4,'Delimiter','\t');
end
fclose(fileID);

% We single out reach 7 for use later on
Reach_7 = Reach{4,1}; % 7VANEM is the reach going direct into Vanemfj. 

no_days = length(Reach{1,1}{1,1}(:,1));
m_start=[1983, 1, 1]; %[1983, 1, 4]
m_stop = m_start + no_days;

%% Unit conversion
sec_to_day = 60*60*24;
liter_to_cube = 1000;
TDP_to_DOP = 0.2;

%% Loading weather

MyLake_weather_only = strcat('C:\OneDrive\PROJECT\O-12111_MARS\4.4_Northern_basin\Vansjø\Model\IO\vansjo_weather_only_',clim_ID,'.txt'); %  this is the weather template pre-made

%% Preparing MyL input
% Thus  combine reaches 5Vaaler and 6Store together
New_Reach = fn_INCA_reach_combination(Reach,2,3,7); % Reach array, ID1, ID2, number of variables in the reach file)., We combine reach 5 and 6 basically. 

inflow = New_Reach{1,1}.*sec_to_day; % m3/sec to m3/day
inflowTemp = New_Reach{1,7};
inflowC = zeros(no_days,1);
inflowS = New_Reach{1,2}.*liter_to_cube; % mg/L to mg/m3
inflowTP = New_Reach{1,5}.*liter_to_cube; % mg/L to mg/m3
inflowDOP = New_Reach{1,3}.*liter_to_cube*TDP_to_DOP; % warning this gives a high humber. is TDP really DIP ? Factor 5
inflowChl = inflowC; % still fixed
inflowDOC = inflowC + 3000; % still fixed
inflowDIC = inflowC; %dummy
inflowO = inflowC; %dummy
%outflowDIC = DICzt(1, :)';
%outflowO = O2zt(1, :)';

INCA_input = [inflow inflowTemp inflowC inflowS inflowTP inflowDOP inflowChl inflowDOC inflowDIC inflowO];

INCA_QC_1 = INCA_input;

store_INCAP_input = tempname;
merge_l_b_inputs(MyLake_weather_only,INCA_input,store_INCAP_input)

%% Writing Vanemfj.
inflow = Reach_7{1,1}.*sec_to_day;
inflowTemp = Reach_7{1,7};
inflowC = zeros(no_days,1);
inflowS = Reach_7{1,2}.*liter_to_cube;
inflowTP = Reach_7{1,5}.*liter_to_cube; %
inflowDOP = Reach_7{1,3}.*liter_to_cube*TDP_to_DOP; % warning this gives a high humber. is TDP really DIP ?
inflowChl = inflowC;
inflowDOC = inflowC + 3000;
inflowDIC = inflowC; %dummy for MyLake TSA
inflowO = inflowC; %dummy for MyLake TSA

INCA_input = [inflow inflowTemp inflowC inflowS inflowTP inflowDOP inflowChl inflowDOC inflowDIC inflowO];

INCA_QC_2 = INCA_input;

vanem_INCAP_input = tempname;
merge_l_b_inputs(MyLake_weather_only,INCA_input,vanem_INCAP_input)

INCA_QC = {INCA_QC_1, INCA_QC_2};

end