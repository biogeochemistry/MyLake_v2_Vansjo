% make text file for MARS all partner meeting according to template:
% Year	month	day	streamflow(mm/day)	air temperature (C)	concentration(ug/l)
% 1973	01      01	2.5                 -5.3                0.025
% tab delimited

%% big_results loop
clear date y m s d test_date day month year

%Big_results(1,X) =
%(:,1) Datum
%(:,2) TP
%(:,3) Chl
%(:,4) PO4
%(:,5) Part
%(:,6) Water temperature
%(:,7) Q_lake

state_var = {...
    'dummy',...
    'TP',...
    'Chl',...
    'PO4',...
    'Ppart',...
    'T'};

scenario = {...
    'NO_Vansjo_Hist_',...
    'NO_Vansjo_Base_',...
    'NO_Vansjo_Hist_M2_',...
    'NO_Vansjo_Hist_M3_',...
    'NO_Vansjo_G4_',...
    'NO_Vansjo_I4_',...
    'NO_Vansjo_G8_',...
    'NO_Vansjo_I8_',...
    'NO_Vansjo_G8Tech_',...
    'NO_Vansjo_G4Cons_',...
    'NO_Vansjo_G8Frag_',...
    'NO_Vansjo_I8Tech_',...
    'NO_Vansjo_I4Cons_',...
    'NO_Vansjo_I8Frag_'};

% loop to nowhere , eventually loop thorugh all vars
for cur_sv = 2:length(state_var);
    for current_run = 1:no_runs
        % CountryCode_Basin_Scenario_Variable.csv is the naming scheme
        % I4 I8, G4, G8, G8Tech, G4Conc
        
        run_ID = strcat(scenario{current_run},state_var(cur_sv)) %  CALIBRATION RUN
        
        % Writing for Task 4.4 common output format
        
        date = big_results{1,current_run}(:,1);
        [y, m, d] = datevec(date,'dd.mm.yy');
        T = big_inputs{1,current_run}(:,3);
        Q = big_results{1,current_run}(:,7) ./ 24 ./ 634;
        % current_state_var
        var_write_out = big_results{1,current_run}(:,cur_sv);
        
        to_write_WP4 = [y, m, d, Q, T, var_write_out];
        
        run_ID_WP4 = strcat(run_ID,'.csv');
        
        cd Postproc_code\MARS
        
        fid=fopen(run_ID_WP4{1},'wt');
        fprintf(fid, '%s \t %s \t %s  \t %s \t %s \t %s', 'Year,',	'Month,'	,'Day,',	'Streamflow(mm/day),',	'Air_Temperature(C),',	'Concentration(ug/l)');
        fprintf(fid,'\n');
        fclose(fid);
        
        dlmwrite(run_ID_WP4{1},to_write_WP4 , '-append');
        cd ..\..
        
        % create aggregated time-series for R scripts
        cd ..\..\..\6.1_Meta_analysis\R_scripts
        % take this array, keep only month (col 2 is >4 and <10) may-sept
        
        
        %M - your data matrix with 4 columns < Year Month Day data > - double type
        [a,~,c] = unique(to_write_WP4(:,1),'rows'); % (:,1) for years, (:,1:2) for months
        accumul_Q = [a, accumarray(c,to_write_WP4(:,4),[],@mean)];
        accumul_T = [a, accumarray(c,to_write_WP4(:,5),[],@mean)];
        accumul_var_writeout = [a, accumarray(c,to_write_WP4(:,6),[],@mean)];
        
        run_ID_to_R = strcat(run_ID,'_yr.csv');
        
        fid=fopen(run_ID_to_R{1},'wt');
        fprintf(fid, '%s \t %s \t %s  \t %s \t %s \t %s', 'Year,','Streamflow(mm/day),',	'Air_Temperature(C),',	'Concentration(ug/l)');
        fprintf(fid,'\n');
        fclose(fid);
        
        to_write_WP6  = [accumul_Q(:,1), accumul_Q(:,2), accumul_T(:,2), accumul_var_writeout(:,2)];
        dlmwrite(run_ID_to_R{1},to_write_WP6 , '-append');
        
        cd ..\..\4.4_Northern_basin\Vansjø\Model
        
    end
end

%% big_inputs(1,X) =
%(:,1)  Global radiation (MJ/(m^2 day))
%(:,2)  Cloud cover (-)
%(:,3)  Air temperature (deg. C, at 2 m height)
%(:,4)  Relative humidity (%, at 2 m height)
%(:,5)  Air pressure (mbar)
%(:,6)  Wind speed (m/s at 10 m height)
%(:,7)  Precipitation (mm/day)
%(:,8)  Inflow volume (m3 day-1)
%(:,9)  Inflow temperature (deg C)
%(:,10) Inflow tracer concentration (-)
%(:,11) Inflow sedimenting tracer (or suspended inorganic matter) concentration (kg m-3)
%(:,12) Inflow total phosphorus (TP) concentration (mg m-3)
%(:,13) Inflow dissolved organic phosphorus (DOP) concentration (mg m-3)
%(:,14) Inflow chlorophyll a concentration (mg m-3)
%(:,15) Inflow DOC concentration (mg m-3)

state_var = ...
    {'Rad',...
    'Cover',...
    'T_air',...
    'Hum',...
    'Patm',...
    'Wind',...
    'Precpi',...
    'Q_stream',...
    'temp_stream',...
    'tracer',...
    'sed_stream',...
    'TP_stream',...
    'DOP_stream',...
    'Chl_stream',...
    'Doc_stream'};

cur_sv = 12; % we want only stream TP

for current_run = 1:no_runs
    % CountryCode_Basin_Scenario_Variable.csv is the naming scheme
    % I4 I8, G4, G8, G8Tech, G4Conc
    
    run_ID = strcat(scenario{current_run},state_var(cur_sv)) %  CALIBRATION RUN
    
    % Writing for Task 4.4 common output format
    
    date = big_results{1,current_run}(:,1);
    [y, m, d] = datevec(date,'dd.mm.yy');
    T = big_inputs{1,current_run}(:,3);
    Q = big_results{1,current_run}(:,7) ./ 24;
    % current_state_var
    var_write_out = big_inputs{1,current_run}(:,cur_sv);
    
    
    to_write_WP4 = [y, m, d, Q, T, var_write_out];
    
    % creating an index vector of summer months
    summer_index = to_write_WP4(:,2)>4 & to_write_WP4(:,2)<9;
    
    run_ID_WP4 = strcat(run_ID,'.csv');
    
    cd Postproc_code\MARS
    
    fid=fopen(run_ID_WP4{1},'wt');
    fprintf(fid, '%s \t %s \t %s  \t %s \t %s \t %s', 'Year,',	'Month,'	,'Day,',	'Streamflow(mm/day),',	'Air_Temperature(C),',	'Concentration(ug/l)');
    fprintf(fid,'\n');
    fclose(fid);
    
    dlmwrite(run_ID_WP4{1},to_write_WP4 , '-append');
    cd ..\..
    
    % create aggregated time-series for R scripts
    cd ..\..\..\6.1_Meta_analysis\R_scripts
    
    %M - your data matrix with 4 columns < Year Month Day data > - double type
    
    to_write_WP4 = to_write_WP4(summer_index,:);
    
    [a,~,c] = unique(to_write_WP4(:,1),'rows'); % (:,1) to retreive years, (:,1:2) for months
    accumul_Q = [a, accumarray(c,to_write_WP4(:,4),[],@mean)];
    accumul_T = [a, accumarray(c,to_write_WP4(:,5),[],@mean)];
    accumul_var_writeout = [a, accumarray(c,to_write_WP4(:,6),[],@mean)];
    
    run_ID_to_R = strcat(run_ID,'_yr.csv');
    
    fid=fopen(run_ID_to_R{1},'wt');
    fprintf(fid, '%s \t %s \t %s  \t %s \t %s \t %s', 'Year,','Streamflow(mm/day),',	'Air_Temperature(C),',	'Concentration(ug/l)');
    fprintf(fid,'\n');
    fclose(fid);
    
    to_write_WP6  = [accumul_Q(:,1), accumul_Q(:,2), accumul_T(:,2), accumul_var_writeout(:,2)];
    
    dlmwrite(run_ID_to_R{1},to_write_WP6 , '-append');
    
    cd ..\..\4.4_Northern_basin\Vansjø\Model
    
end
