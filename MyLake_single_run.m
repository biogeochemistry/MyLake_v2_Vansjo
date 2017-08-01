% for i=1:1000
tic
disp('Started at:')
disp(datetime('now'));
[lake_params, sediment_params] = load_params();

run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

no_runs = 1; % 26/7/2016 ... did not find no_run so I added it again
is_metrics = true; % print metrics in the end

big_results = cell(1,no_runs);  % collects the results
big_inputs = cell(1,no_runs);   % collects the inputs
is_save_results = false;


% Initial calibration for P:
x = [0.0643, 0.1368, 1.4647, 0.5293, 0.0570, 0.1000, 1.3339, 0.8793, 0.0590, 101.1213, 68.1137]; % ecomac-2 final result, RMSD 44.0151
% x = [0.3050, 0.2560, 1.1311, 0.4697, 0.1633, 0.1002, 1.4969, 1.5477, 0.5041, 6.754e+5, 94.3578]; % 27.07.2017 RMSD 56


lake_params{52 -5} = x(1); % 9     settling velocity for Chl1 a (m day-1)
lake_params{54 -5} = x(2); % 11    loss rate (1/day) at 20 deg C
lake_params{55 -5} = x(3); % 12    specific growth rate (1/day) at 20 deg C
lake_params{58 -5} = x(4); % 15    Half saturation growth P level (mg/m3)
lake_params{61 -5} = x(5); % 18    Settling velocity for Chl2 a (m day-1)
lake_params{62 -5} = x(6);  % 19    Loss rate (1/day) at 20 deg C
lake_params{63 -5} = x(7);  % 20    Specific growth rate (1/day) at 20 deg C
lake_params{64 -5} = x(8);  % 21    Half saturation growth P level (mg/m3)
lake_params{51 -5} = x(9);  % % 8  settling velocity for S (m day-1)
sediment_params{23} = x(10);  % 38 R16 sorption of P on Fe k
sediment_params{35} = x(11);  %    accel

sediment_params{35} = 50;  %    accel
lake_params{51 -5} = 0.5;  %    % 8   settling velocity for S (m day-1)
lake_params{61 -5} = 0.2; % 18    Settling velocity for Chl2 a (m day-1)
lake_params{52 -5} = 0.2; % 9     settling velocity for Chl1 a (m day-1)
sediment_params{23} = 1e1;  % 38 R16 sorption of P on Fe k

% % Temperature calibration results:
lake_params{2} =   44.2401e-003;
lake_params{5} = 500.0000e-003;
% lake_params{5} = 0.3;


% % Scaling factors
lake_params{21 -5} = 1; % 16    scaling factor for inflow volume (-)
lake_params{22 -5} = 1; % 17    adjusting delta for inflow temperature (-)
lake_params{23 -5} = 1; % 18    scaling factor for inflow concentration of C (-)
lake_params{24 -5} = 1; % 19    scaling factor for inflow concentration of S (-)
lake_params{25 -5} = 1; % 20    scaling factor for inflow concentration of total P (-)
lake_params{26 -5} = 1; % 21    scaling factor for inflow concentration of diss. organic P (-)
lake_params{27 -5} = 1; % 22    scaling factor for inflow concentration of Chl a (-)
lake_params{28 -5} = 1; % 23    scaling factor for inflow concentration of DOC  (-)
lake_params{29 -5} = 1; % 24    scaling factor for inflow concentration of POC  (-)
lake_params{30 -5} = 1; % 25    Scaling factor for inflow concentration of O2 (-)
lake_params{31 -5} = 1; % 26    Scaling factor for inflow concentration of DOC  (-)
lake_params{32 -5} = 1; % 27    Scaling factor for inflow concentration of NO3 (-)
lake_params{33 -5} = 1; % 28    Scaling factor for inflow concentration of NH4 (-)
lake_params{34 -5} = 1; % 29    Scaling factor for inflow concentration of SO4 (-)
lake_params{35 -5} = 1; % 30    Scaling factor for inflow concentration of Fe2 (-)
lake_params{36 -5} = 1; % 31    Scaling factor for inflow concentration of Ca2 (-)
lake_params{37 -5} = 1; % 32    Scaling factor for inflow concentration of pH (-)
lake_params{38 -5} = 1; % 33    Scaling factor for inflow concentration of CH4 (-)
lake_params{39 -5} = 1; % 34    Scaling factor for inflow concentration of Fe3 (-)
lake_params{40 -5} = 1; % 35    Scaling factor for inflow concentration of Al3 (-)
lake_params{41 -5} = 1; % 36    Scaling factor for inflow concentration of SiO4 (-)
lake_params{42 -5} = 1; % 37    Scaling factor for inflow concentration of SiO2 (-)
lake_params{43 -5} = 1; % 38    Scaling factor for inflow concentration of diatom (-)


% O&P calibration:
x = [316.9331e-003, 212.4112e-003, 1.3954e+000, 964.7451e-003, 500.0000e-003, 112.2803e-003, 1.4675e+000, 1.4823e+000, 757.5944e-003, 1.0000e+000, 1.0000e+000, 901.6600e-003, 8.1686e+000, 1.2091e+000, 221.2731e-003, 2.7304e+000, 0.0000e+000];

x = [500.0000e-003,   101.0838e-003,     1.0000e+000,     1.6996e+000,   188.6713e-003,   100.0000e-003,     1.4400e+000,     1.6311e+000,    10.0000e-003,     1.0000e+000,     2.0000e+000,     2.7065e+000,    4*62.1176e-003,    85.8779e-003]; % RMSD 137.82

lake_params{52 -5} = x(1); % 9     settling velocity for Chl1 a (m day-1)
lake_params{54 -5} = x(2); % 11    loss rate (1/day) at 20 deg C
lake_params{55 -5} = x(3); % 12    specific growth rate (1/day) at 20 deg C
lake_params{58 -5} = x(4); % 15    Half saturation growth P level (mg/m3)
lake_params{61 -5} = x(5); % 18    Settling velocity for Chl2 a (m day-1)
lake_params{62 -5} = x(6);  % 19    Loss rate (1/day) at 20 deg C
lake_params{63 -5} = x(7);  % 20    Specific growth rate (1/day) at 20 deg C
lake_params{64 -5} = x(8);  % 21    Half saturation growth P level (mg/m3)
lake_params{51 -5} = x(9);  % % 8  settling velocity for S (m day-1)
sediment_params{23} = x(10);  % 38 R16 sorption of P on Fe k
sediment_params{35} = x(11);  %    accel
% lake_params{21 -5} = x(12); % 16    scaling factor for inflow volume (-)
% lake_params{25 -5} = x(13); % 20    scaling factor for inflow concentration of total P (-)
% lake_params{26 -5} = x(14); % 21    scaling factor for inflow concentration of diss. organic P (-)
lake_params{28 -5} = x(12); % 23    scaling factor for inflow concentration of DOC  (-)
lake_params{29 -5} = x(13); % 23    scaling factor for inflow concentration of POC  (-)
lake_params{39 -5} = x(14); % 34    Scaling factor for inflow concentration of Fe3 (-)


% % Trials:
lake_params{24 -5} = 100; % 23    scaling factor for inflow concentration of POC  (-)
lake_params{51 -5} = 0.01;  % % 8  settling velocity for S (m day-1)
% lake_params{15 -5} = 2e-4; % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{59 -5} = 1e-5; % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{17 -5} = 0.015; % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{60 -5} = 0.015; % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{35} = 50;  %    accel

lake_params{61 -5} = 0.01; % 18    Settling velocity for Chl2 a (m day-1)
% lake_params{62 -5} = 0.1;  % 19    Loss rate (1/day) at 20 deg C
% lake_params{63 -5} = 2;  % 20    Specific growth rate (1/day) at 20 deg C
% lake_params{64 -5} = 6;  % 21    Half saturation growth P level (mg/m3)

% parfor
for current_run = 1:no_runs
    if current_run == 1;
        run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
        clim_ID = run_ID
        if use_INCA == 1
            m_start=[2000, 1, 1]; %
            m_stop=[2012, 12, 31]; %
        else
            m_start=[2004, 1, 1]; %
            m_stop=[2010, 12, 31]; %
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
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, is_save_results); % runs the model and outputs obs and sim
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

% if is_metrics == true

%     zinx=find(MyLake_results.basin1.z<4);
%     TP_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:)+MyLake_results.basin1.concentrations.PP(zinx,:) + MyLake_results.basin1.concentrations.Chl(zinx,:)+MyLake_results.basin1.concentrations.C(zinx,:)+MyLake_results.basin1.concentrations.DOP(zinx,:))', 2);
%     Chl_mod = mean((MyLake_results.basin1.concentrations.Chl(zinx,:)+MyLake_results.basin1.concentrations.C(zinx,:))', 2);
%     P_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:))', 2);
%     PP_mod = mean((MyLake_results.basin1.concentrations.PP(zinx,:))', 2);

%     load 'obs/store_obs/TOTP.dat' % measured
%     % load 'obs/store_obs/Cha.dat' % measured
%     load 'obs/store_obs/Cha_aquaM_march_2017.dat' % measured
%     load 'obs/store_obs/PO4.dat' % measured
%     load 'obs/store_obs/Part.dat' % measured


%     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, TOTP(:,1)));
%     r_TOTP = RMSE(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


%     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha_aquaM_march_2017(:,1)));
%     r_Chl = RMSE(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));


%     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
%     r_PO4 = RMSE(P_mod(loc_sim, 1), PO4(loc_obs, 2));


%     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
%     r_PP = RMSE(PP_mod(loc_sim, 1), Part(loc_obs, 2));

%     disp('RMSD P:')
%     disp(sum([r_TOTP, r_Chl, r_PO4, r_PP]))
% end


toc




% end
