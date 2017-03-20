tic
[lake_params, sediment_params] = load_params();

run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 1; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

no_runs = 1 % 26/7/2016 ... did not find no_run so I added it again

big_results = cell(1,no_runs);  % collects the results
big_inputs = cell(1,no_runs);   % collects the inputs


% lake_params{36-7} =2.5196*1e3; %  2500    Half saturation parameter for Langmuir isotherm
% lake_params{37-7} =7.3827*1e3; %   8000    Scaling parameter for Langmuir isotherm !!!!!!!!!!!!
% lake_params{39-7} =0.00002*1e3; %   0.02    settling velocity for Chl1 a (m day-1)
% lake_params{41-7} =0.0001*1e3; %   0.2    loss rate (1/day) at 20 deg C
% lake_params{42-7} =0.0020*1e3; %   1.5    specific growth rate (1/day) at 20 deg C
% lake_params{45-7} =0.0002*1e3; %   0.2    Half saturation growth P level (mg/m3)
% lake_params{48-7} =0.00002*1e3; %   0.02    Settling velocity for Chl2 a (m day-1)
% lake_params{49-7} =0.0002*1e3;  %    0.2    Loss rate (1/day) at 20 deg C
% lake_params{50-7} =0.0008*1e3;  %   1.5    Specific growth rate (1/day) at 20 deg C
% lake_params{51-7} =0.0003*1e3;  %   0.2    Half saturation growth P level (mg/m3)
% lake_params{68-7} =0.0023*1e3;  %   2];    Q10 for reactions of respiration





% parfor
for current_run = 1:no_runs
    if current_run == 1;
        run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
        clim_ID = run_ID
        if use_INCA == 1
            m_start=[2005, 1, 1]; %
            m_stop=[2005, 12, 31]; %
        else
            m_start=[2000, 1, 1]; %
            m_stop=[2012, 12, 31]; %
        end

    elseif current_run == 2;
        run_ID = 'Vansjo_Hist_M1' ; % NOT REAL WEATHER FOR THIS ONE !
        clim_ID = run_ID
        m_start=[1995, 1, 1]; %
        m_stop=[2070, 12, 31]; %

    elseif current_run == 3;
        run_ID = 'Vansjo_Hist_M2' ; %
        clim_ID = 'Vansjo_Hist_M0' ;
        m_start=[1984, 1, 1]; %
        m_stop=[2013, 12, 31]; %

    elseif current_run == 4;
        run_ID = 'Vansjo_Hist_M3' ; %
        clim_ID = 'Vansjo_Hist_M0' ;
        m_start=[1984, 1, 1]; %
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


    % try
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID); % runs the model and outputs obs and sim
        big_results{current_run} = {MyLake_results, Sediment_results};
    % catch ME
    %     fprintf('Process crashed: %s\n', num2str(current_run))
    %     fprintf('\tID: %s\n', ME.identifier)
    %     fprintf('\tMessage: %s\n', ME.message)
    %     fprintf('\tStack::\n')
    %     for k=1:length(ME.stack)
    %         disp(ME.stack(k))
    %     end
    % end

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
