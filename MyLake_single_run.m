for i=1:1000
tic
disp('Started at:')
disp(datetime('now'));
[lake_params, sediment_params] = load_params();

run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

no_runs = 1; % 26/7/2016 ... did not find no_run so I added it again

big_results = cell(1,no_runs);  % collects the results
big_inputs = cell(1,no_runs);   % collects the inputs

% x = [0.00365275111636845, 0.444888225060776, 6.21096051969837, 1.98304729966760, 0.00423770963447693, 0.585402600080357, 4.17074391909280, 1.84803487767737, 0.0420205265067727, 84.9948291324757, 30];

% x = [0.0957, 0.8462, 12.2918, 0.4898, 0.1904, 1.3158, 9.7200, 1.4780, 0.02, 10.0220, 40.9828];
% x = [0.4958, 0.1998, 1.1923, 1.4176, 0.1125, 0.2737, 1.3935, 0.8685, 0.6889, 259.7042, 40.3020];
% x = [0.4958, 0.1054, 1.3664, 1.7035, 0.0698, 0.1546, 1.3697, 0.9468, 0.02, 100.7042, 24.7989];

% x = [245.0402e-003, 193.8201e-003, 1.3948e+000, 700.8806e-003, 76.1493e-003, 121.8722e-003, 1.3160e+000, 363.7310e-003, 77.1076e-003, 443.2972e+000, 34.9491e+000];

% new Lang test
% x = [0.4958, 0.1054, 1.3664, 1.7035, 0.0698, 0.1546, 1.3697, 0.9468, 0.2, 1e-3, 24.7989];

% x = [0.0643, 0.1368, 1.1846, 0.2516, 0.0570, 0.1000, 1.3924, 0.5457, 0.0590, 101.1213, 67.8637]; x'; % ecomac-2, best at the moment
% x = [0.0643, 0.1368, 1.4647, 0.5293, 0.0570, 0.1000, 1.3339, 0.8793, 0.0590, 101.1213, 67.8637] % ecomac-2, rmsd 44.01
x = [0.0643, 0.1368, 1.4647, 0.5293, 0.0570, 0.1000, 1.3339, 0.8793, 0.0590, 101.1213, 68.1137]; % ecomac-2 final result, RMSD 44.0151
% x = [0.0643, 0.1800, 1.3659, 0.5293, 0.0570, 0.1000, 1.3339, 0.5457, 0.0590, 635.5627, 46.1497]; x'; % ecomac-2 rmsd 44.02
% x = [173.8356e-003, 105.3827e-003, 1.1923e+000, 1.7035e+000, 69.8043e-003, 154.5866e-003, 1.3898e+000, 946.7929e-003, 366.6823e-003, 253.7042e+000, 19.3020e+000]; % ecomac-1, final result, 32.02 final result

lake_params{40 -7} = x(1); % 9     settling velocity for Chl1 a (m day-1)
lake_params{42 -7} = x(2); % 11    loss rate (1/day) at 20 deg C
lake_params{43 -7} = x(3); % 12    specific growth rate (1/day) at 20 deg C
lake_params{46 -7} = x(4); % 15    Half saturation growth P level (mg/m3)
lake_params{49 -7} = x(5); % 18    Settling velocity for Chl2 a (m day-1)
lake_params{50 -7} = x(6);  % 19    Loss rate (1/day) at 20 deg C
lake_params{51 -7} = x(7);  % 20    Specific growth rate (1/day) at 20 deg C
lake_params{52 -7} = x(8);  % 21    Half saturation growth P level (mg/m3)
lake_params{39 -7} = x(9);  %       settling velocity for S (m day-1)
sediment_params{22} = x(10);  % 38 R16 sorption of P on Fe k
sediment_params{34} = x(11);  %    accel

sediment_params{34} = 1;  %    accel
lake_params{39 -7} = 0.5;  %       settling velocity for S (m day-1)
lake_params{49 -7} = 0.1; % 18    Settling velocity for Chl2 a (m day-1)
lake_params{40 -7} = 0.1; % 9     settling velocity for Chl1 a (m day-1)

% parfor
for current_run = 1:no_runs
    if current_run == 1;
        run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
        clim_ID = run_ID
        if use_INCA == 1
            m_start=[2000, 1, 1]; %
            m_stop=[2012, 12, 31]; %
        else
            m_start=[2000, 1, 1]; %
            m_stop=[2013, 12, 31]; %
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


disp('Saving results...')
save('IO/MyLakeResults.mat', 'MyLake_results', 'Sediment_results')
disp('Finished at:')
disp(datetime('now'));
toc
end
