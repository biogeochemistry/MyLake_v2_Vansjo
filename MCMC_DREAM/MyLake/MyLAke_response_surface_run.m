tic
MyLake_setPriors_for_MCMC % prepares parameter files
clear big_result
m_start=[1984, 5, 13]; % 1973 1 1 is the earliest possible given the provided input files8
m_stop=[1986, 5, 31]; % 2012, 12, 31 is the latest possible given the provided input files  

run_INCA = 0; % 1- MyLake will run INCA, 0- No run

big_result = zeros(2,length(list_input_file),(datenum(m_stop)-datenum(m_start)));

parfor current_run = 1:length(list_input_file)

            current_filename = list_input_file{current_run};
            use_INCA=tempname; % we trick the MyLake Application to get the inputfilename as the use_INCA parameter ... 
            copyfile(current_filename,use_INCA);
            % run model and pass inputfile as an argument ... 
            [TP_obs,TP_mod,chl_obs,chl_mod] = fn_MyL_application(m_start, m_stop, K_values_sediments, K_values_lake, use_INCA, run_INCA); % runs the model and outputs obs and sim
            
            % accumulating results
            big_result(1,current_run,:) = TP_mod(:,2); 
            big_result(2,current_run,:) = chl_mod(:,2);
            % clean up
            delete(use_INCA)
end

% reformmating into 3D array

 toc