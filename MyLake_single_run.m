tic
load_params;

run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 1; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

no_runs = 1 % 26/7/2016 ... did not find no_run so I added it again

% for calibratoin against real data, set all inca to 0 and no_runs to 1.
% then, disable the "big_result" array down.

big_results = cell(1,no_runs);  % collects the results
big_inputs = cell(1,no_runs);   % collects the inputs

% Scenario = {...
% 'NO_Vansjo_Hist_',...
% 'NO_Vansjo_Base_'...
% 'NO_Vansjo_Hist_M2_'...
% 'NO_Vansjo_Hist_M3_'...
% 'NO_Vansjo_G4_',...
% 'NO_Vansjo_I4_',...
% 'NO_Vansjo_G8_',...
% 'NO_Vansjo_I8_',...
% 'NO_Vansjo_G8Tech_',...
% 'NO_Vansjo_G4Cons_',...
% 'NO_Vansjo_G8Frag_',...
% 'NO_Vansjo_I8Tech_',...
% 'NO_Vansjo_I4Cons_',...
% 'NO_Vansjo_I8Frag_',};

% parfor
for current_run = 1:no_runs
    if current_run == 1;
        run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
        clim_ID = run_ID
        if use_INCA == 1
            m_start=[2004, 1, 1]; % for scenario runs
            m_stop=[2009, 12, 31]; % for scenario runs
        else
            m_start=[2004, 1, 1]; %
            m_stop=[2009, 12, 31]; %
        end

    elseif current_run == 2;
        run_ID = 'Vansjo_Hist_M1' ; % NOT REAL WEATHER FOR THIS ONE !
        clim_ID = run_ID
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 3;
        run_ID = 'Vansjo_Hist_M2' ; %
        clim_ID = 'Vansjo_Hist_M0' ;
        m_start=[1983, 1, 1]; %
        m_stop=[2013, 12, 31]; %

    elseif current_run == 4;
        run_ID = 'Vansjo_Hist_M3' ; %
        clim_ID = 'Vansjo_Hist_M0' ;
        m_start=[1983, 1, 1]; %
        m_stop=[2013, 12, 31]; %

    elseif current_run == 5;
        run_ID = 'Vansjo_RCP4_GFDL_M0' ; %
        clim_ID = run_ID
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 6;
        run_ID = 'Vansjo_RCP4_IPSL_M0' ; %
        clim_ID = run_ID
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 7;
        run_ID = 'Vansjo_RCP8_GFDL_M0' ; %
        clim_ID = run_ID
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 8;
        run_ID = 'Vansjo_RCP8_IPSL_M0' ; %
        clim_ID = run_ID
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 9;
        run_ID = 'Vansjo_RCP8_GFDL_M4' ;
        clim_ID = 'Vansjo_RCP8_GFDL_M0' ; %
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 10;
        run_ID = 'Vansjo_RCP4_GFDL_M5' ;
        clim_ID = 'Vansjo_RCP4_GFDL_M0' ; %
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 11;
        run_ID = 'Vansjo_RCP8_GFDL_M6' ;
        clim_ID = 'Vansjo_RCP8_GFDL_M0' ; %
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 12;
        run_ID = 'Vansjo_RCP8_IPSL_M4' ;
        clim_ID = 'Vansjo_RCP8_IPSL_M0' ; %
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 13;
        run_ID = 'Vansjo_RCP4_IPSL_M5' ;
        clim_ID = 'Vansjo_RCP4_IPSL_M0' ; %
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    else
        run_ID = 'Vansjo_RCP8_IPSL_M6' ;
        clim_ID = 'Vansjo_RCP8_IPSL_M0' ; %
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    end


    try
        [TP_obs,TP_mod, TP_date,chl_obs,chl_mod, Chl_date, PO4_obs, PO4_mod, PO4_date, Part_obs, Part_mod, Part_date, MyLake_results, Sediment_results, input_all]  = fn_MyL_application(m_start, m_stop, K_values_sediment, K_values_lake, use_INCA, run_INCA, run_ID, clim_ID); % runs the model and outputs obs and sim
        big_results{current_run} = {TP_obs,TP_mod, TP_date,chl_obs,chl_mod, Chl_date, PO4_obs, PO4_mod, PO4_date, Part_obs, Part_mod, Part_date, MyLake_results, Sediment_results, input_all}
    catch ME
        fprintf('Process crashed: %s\n', num2str(current_run))
        fprintf('\tID: %s\n', ME.identifier)
        fprintf('\tMessage: %s\n', ME.message)
        fprintf('\tStack::')
        disp(ME.stack(1))
    end

end

% cd .. ; cd .. ;

% big_results(2,:) = {'Hist_M0','Hist_M1','Hist_M2','Hist_M3','RCP4_GFDL_M0','RCP4_ISPL_M0','RCP8_GFDL_M0','RCP8_ISPL_M0','RCP8_GFDL_M4','RCP4_GFDL_M5','RCP8_GFDL_M6','RCP8_IPSL_M4','RCP4_IPSL_M5','RCP8_IPSL_M6',};
% big_inputs(2,:) = {'Hist_M0','Hist_M1','Hist_M2','Hist_M3','RCP4_GFDL_M0','RCP4_ISPL_M0','RCP8_GFDL_M0','RCP8_ISPL_M0','RCP8_GFDL_M4','RCP4_GFDL_M5','RCP8_GFDL_M6','RCP8_IPSL_M4','RCP4_IPSL_M5','RCP8_IPSL_M6',};

% dev1 = TP_mod(:)-TP_obs(:);
% dev1(isnan(dev1))=[]; %removing NaNs
% SS(1) = nansum((TP_mod(:)-TP_obs(:)).^2);
% Nobs(1) = numel(dev1);
%
% dev2 = chl_mod(:)-chl_obs(:);
% dev2(isnan(dev2))=[];
% SS(2) = nansum((chl_mod(:)-chl_obs(:)).^2);
% Nobs(2) = numel(dev2);
%
% clear dev1 dev2 Nobs % note: cannot clear within parfor

toc
