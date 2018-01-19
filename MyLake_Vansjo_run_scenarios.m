scen = {'T_only_full_scen_base_historical_20y', 'T_only_RCP4_GFDL', 'T_only_RCP4_IPSL', 'T_only_RCP8_GFDL', 'T_only_RCP8_IPSL', 'T_only_RCP45_NorESM', 'T_only_RCP60_NorESM'}

parfor s = 1:5

    name_of_scenario = strcat('IO/airT_Scenarios/', scen{s}, '.txt')
    file_name = strcat('IO/airT_Scenarios/96ts_', scen{s}, '.mat')



    tic
    disp('Started at:')
    disp(datetime('now'));

    is_metrics = true; % print metrics in the end

    m_start=[1995, 1, 1]; %
    m_stop=[2050, 12, 31]; %
    % big_results = cell(1,no_runs);  % collects the results
    % big_inputs = cell(1,no_runs);   % collects the inputs
    save_initial_conditions = false; % save final concentrations as initial for the next run


    [lake_params, sediment_params] = load_params();
    % name_of_scenario = 'IO/airT_Scenarios/T_only_RCP4_IPSL.txt'
    % file_name = 'IO/airT_Scenarios/T_only_RCP4_IPSL.mat'


    % inaccurate but faster:
    sediment_params{73} = 96;
    sediment_params{74} = 0; % pH algo disabled;
    % sediment_params{72} = 0; % effective depth test

    % ======================================================================
    % file_name = 'IO/test.mat'
    % % chl calibration;
    % lake_params{47} = 1.1134e-01; 1.8772e-01; % 3.4430e-01; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
    % lake_params{49} =    1.1534e-01; 1.9668e-01; % 1.3488e-01; % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
    % lake_params{50} =    1.4769e+00; 1.4858e+00; % 1.4902e+00; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
    % lake_params{53} =    5.3297e-01; 1.7330e+00; % 9.5346e-01; % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
    % lake_params{56} =    3.6577e-01; 8.0391e-02; % 2.3220e-01; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
    % lake_params{57} =    1.3731e-01; 1.6195e-01; % 2.1678e-01; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
    % lake_params{58} =    1.4592e+00; 1.4877e+00; % 1.4685e+00; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
    % lake_params{59} =    8.0115e-01; 1.0889e+00; % 6.5863e-01; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
    % lake_params{46} =    1.3786e-01; 1.1327e-01; % 2.8629e-02; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
    % lake_params{10} =    8.7271e-05; 8.6703e-05; % 6.8345e-05;  % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    % lake_params{54} =    4.6981e-05; 5.0520e-05; % 5.0846e-05;  % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    % lake_params{12} =    4.4221e-02; 4.5000e-02; % 2.0724e-02;  % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
    % lake_params{55} =    3.8198e-02; 4.5000e-02; % 1.9385e-02;  % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
    % sediment_params{52} =    5.2754e+00; 2.0000e+00;  % 65.1237e+000   %    accel
    % lake_params{24} = 5.0100e-01;  % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)

    % % trial
    % lake_params{24} = 1.0;  % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
    % lake_params{20} = 1.4;  % 390.1162e-003   % 24    scaling factor for inflow concentration of TP (-)
    % lake_params{46} = 0.08; % 2.8629e-02; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
    % lake_params{72} = 35; %   30  % 72           % depth below which the lake is affected by sediments, [m], if -1 (experimental) , then sediments below pycnocline

    % NIVA latest RMSD only result (++++Best so far++++):
    x = [0.0865681634702761; 0.161705696258267; 1.23114894752963; 1.66632346919548; 0.217316646052962; 0.184295158144136; 1.50000000000000; 1.44746872184578; 0.0609313437550716; 2.26727837370864e-05; 3.06302033352011e-05; 0.0450000000000000; 0.0321267784864950; 19.3181501625625; 1.15335807117479; 0.718058521202605];
    % ======================================================================
    % file_name = 'IO/chl_rmsd.mat'
    % chl calibration RMSD
    lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
    lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
    lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
    lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
    lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
    lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
    lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
    lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
    lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
    lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
    lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
    sediment_params{52} = x(14); % 65.1237e+000   %    accel
    lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
    lake_params{20} = x(16); % 1   % 20    scaling factor for inflow concentration of TP (-)
    % ======================================================================

    % ======================================================================
    % ecomac-1 chl nrmsd ~= 5.
    % file_name = 'IO/chl_nrmsd.mat'
    % x = [0.0579239774149386; 0.116422281177389; 1.38484320007442; 1.53201985840633; 0.160917008960281; 0.178095709943439; 1.26619686863682; 1.58643475878890; 0.146673515447364; 6.78211244339783e-05; 2.17328251518974e-05; 0.0445219170323908; 0.0267633321331395; 11.9400654956391; 1.68021340701733; 1.63593615750069];
    % % chl calibration NRMSD
    % lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
    % lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
    % lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
    % lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
    % lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
    % lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
    % lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
    % lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
    % lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
    % lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    % lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    % lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
    % lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
    % sediment_params{52} = x(14); % 65.1237e+000   %    accel
    % lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
    % lake_params{20} = x(16); % 1   % 20    scaling factor for inflow concentration of TP (-)
    % ======================================================================


    % ======================================================================
    % ecomac-2 sediment cores nrmsd = 41.
    % x = [10.6160463423489; 0.573135424263471; 0.0979950894572524; 0.0621489725692570; 0.0158366788053187; 0.0614849072158739; 0.116104561853594; 83.5585559947998; 1.57839693972576; 2.76943702988488; 6.15906670539647; 1.20706687750177; 35.9606317972236; 63.5661388861377; 0.999195026938994; 10.1533888812312; 2.22460465116220; 7.15678095558815; 0.720516440154948; 1.45372783090552];

    % file_name = 'IO/test_0.mat'
    % % cores calibration
    % sediment_params{1} = x(1); % 65.1237e+000   %    accel
    % sediment_params{1} = x(2);  %   'k_Chl',                 %        % 1
    % sediment_params{2} = x(3);  %  'k_POP',                 %        % 1
    % sediment_params{3} = x(4);  % 'k_POC',                  %        % 0.01
    % sediment_params{4} = x(5);  %  'k_DOP',                 %        % 1
    % sediment_params{5} = x(6);  % 'k_DOC',                  %        % 1
    % sediment_params{23} = x(7);  %     'k_pdesorb_a',         %
    % sediment_params{24} = x(8);  %     'k_pdesorb_b',         %
    % sediment_params{54} = x(9);  %     'k_pdesorb_c',         %
    % lake_params{75} = x(10);%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo
    % lake_params{22} = x(11);%    scaling factor for inflow concentration of Chl a (-)
    % lake_params{25} = x(12);%    Scaling factor for inflow concentration of O2 (-)
    % lake_params{27} = x(13);%    Scaling factor for inflow concentration of NO3 (-)
    % lake_params{34} = x(14);%    Scaling factor for inflow concentration of Fe3 (-)
    % lake_params{35} = x(15);%    Scaling factor for inflow concentration of Al3 (-)
    % lake_params{37} = x(16);%    Scaling factor for inflow concentration of CaCO3 (-)
    % lake_params{31} = x(17);%    k_apa_pre
    % lake_params{32} = x(18);%    k_apa_pre
    % lake_params{40} = x(19);%    k_viv_pre
    % lake_params{41} = x(20);%    k_viv_pre
    % ======================================================================



    % ======================================================================
    % Niva res. NIVA: calibration of sed. with ts=24 and custom weights (3x,2x) + nrmsd + pH 8
    % res ~= 132
    % file_name = 'IO/niva_pH_8_RMSD_chl.mat'

    % calib_res = [0.0865681634702761; 0.161705696258267; 1.23114894752963; 1.66632346919548; 0.217316646052962; 0.184295158144136; 1.50000000000000; 1.44746872184578; 0.0609313437550716; 2.26727837370864e-05; 3.06302033352011e-05; 0.0450000000000000; 0.0321267784864950; 19.3181501625625; 1.15335807117479; 0.718058521202605];

    % % Chl better predictions:
    % lake_params{47} = calib_res(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
    % lake_params{49} = calib_res(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
    % lake_params{50} = calib_res(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
    % lake_params{53} = calib_res(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
    % lake_params{56} = calib_res(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
    % lake_params{57} = calib_res(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
    % lake_params{58} = calib_res(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
    % lake_params{59} = calib_res(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
    % lake_params{46} = calib_res(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
    % lake_params{10} = calib_res(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    % lake_params{54} = calib_res(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    % lake_params{12} = calib_res(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
    % lake_params{55} = calib_res(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
    % sediment_params{52} = calib_res(14); % 65.1237e+000   %    accel
    % lake_params{24} = calib_res(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
    % lake_params{20} = calib_res(16); % 1   % 20    scaling factor for inflow concentration of TP (-)

    % x = [62.5462980835090; 0.761670173737531; 0.0991459596796610; 0.0611300033799665; 0.0555712909906529; 0.00440011130189627; 1.89391344140936; 33.5989496699941; 4.36261755250612; 74.3796267147851; 1.49837337530373; 0.492764738313926; 69.4891706612324; 2.35993627709247; 11.8785780178841; 3.93463786921662; 7.20726205525980; 8.43667914295516; 2.51650476072176; 3.47040978594033];

    % % new added for cores
    % sediment_params{1} = x(1); % 65.1237e+000   %    accel
    % sediment_params{1} = x(2);  %   'k_Chl',                 %        % 1
    % sediment_params{2} = x(3);  %  'k_POP',                 %        % 1
    % sediment_params{3} = x(4);  % 'k_POC',                  %        % 0.01
    % sediment_params{4} = x(5);  %  'k_DOP',                 %        % 1
    % sediment_params{5} = x(6);  % 'k_DOC',                  %        % 1
    % sediment_params{23} = x(7);  %     'k_pdesorb_a',         %
    % sediment_params{24} = x(8);  %     'k_pdesorb_b',         %
    % sediment_params{54} = x(9);  %     'k_pdesorb_c',         %

    % % SO4 boundary
    % lake_params{75} = x(10);%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo

    % % for cores too (scaling unknown inputs):
    % lake_params{22} = x(11);%    scaling factor for inflow concentration of Chl a (-)
    % lake_params{25} = x(12);%    Scaling factor for inflow concentration of O2 (-)
    % lake_params{27} = x(13);%    Scaling factor for inflow concentration of NO3 (-)
    % lake_params{34} = x(14);%    Scaling factor for inflow concentration of Fe3 (-)
    % lake_params{35} = x(15);%    Scaling factor for inflow concentration of Al3 (-)
    % lake_params{37} = x(16);%    Scaling factor for inflow concentration of CaCO3 (-)

    % % P minerals:
    % lake_params{31} = x(17);%    k_apa_pre
    % lake_params{32} = x(18);%    k_apa_pre
    % lake_params{40} = x(19);%    k_viv_pre
    % lake_params{41} = x(20);%    k_viv_pre
    % ======================================================================

    % Niva res. NIVA: calibration of sed. with ts=24 and custom weights (1x) + nrmsd + pH 8, no POP
    % res ~= 132
    % file_name = 'IO/niva_pH_8_NRMSD_chl_1x_weights.mat'
    % x = [9.9259e+01, 1.0000e-02, 1.0000e-03, 7.4466e-02, 8.8164e-02, 1.5672e-03, 1.9609e+00, 4.0020e+01, 1.4153e-01, 4.7362e+01, 1.9585e+01, 0, 0, 6.1781e-01, 1.2500e+00, 2.5000e+00, 2.4012e+00, 1.2781e+00, 1.6830e+00, 9.2443e+00];

    % sediment_params{1} = x(1); % 65.1237e+000   %    accel
    % sediment_params{1} = x(2);  %   'k_Chl',                 %        % 1
    % sediment_params{2} = x(3);  %  'k_POP',                 %        % 1
    % sediment_params{3} = x(4);  % 'k_POC',                  %        % 0.01
    % sediment_params{4} = x(5);  %  'k_DOP',                 %        % 1
    % sediment_params{5} = x(6);  % 'k_DOC',                  %        % 1
    % sediment_params{23} = x(7);  %     'k_pdesorb_a',         %
    % sediment_params{24} = x(8);  %     'k_pdesorb_b',         %
    % sediment_params{54} = x(9);  %     'k_pdesorb_c',         %

    % % SO4 boundary
    % lake_params{75} = x(10);%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo

    % % for cores too (scaling unknown inputs):
    % lake_params{22} = x(11);%    scaling factor for inflow concentration of Chl a (-)
    % lake_params{25} = x(12);%    Scaling factor for inflow concentration of O2 (-)
    % lake_params{27} = x(13);%    Scaling factor for inflow concentration of NO3 (-)
    % lake_params{34} = x(14);%    Scaling factor for inflow concentration of Fe3 (-)
    % lake_params{35} = x(15);%    Scaling factor for inflow concentration of Al3 (-)
    % lake_params{37} = x(16);%    Scaling factor for inflow concentration of CaCO3 (-)

    % % P minerals:
    % lake_params{31} = x(17);%    k_apa_pre
    % lake_params{32} = x(18);%    k_apa_pre
    % lake_params{40} = x(19);%    k_viv_pre
    % lake_params{41} = x(20);%    k_viv_pre
    % % ======================================================================


    % try
    run_ID = 0;
    clim_ID = 0;
    run_INCA = 0; % 1- MyLake will run INCA, 0- No run
    use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

    [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, name_of_scenario, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim





    disp('Saving results...')
    parsave(file_name, MyLake_results, Sediment_results)
    disp('Finished at:')
    disp(datetime('now'));

    % if is_metrics == true

    %     load('Postproc_code/Vansjo/VAN1_data_2017_02_28_10_55.mat')

    %     depths = [5;10;15;20;25;30;35;40];
    %     rmsd_O2 = 0;


    %     for i=1:size(depths,1)
    %         d = depths(i);
    %         zinx=find(MyLake_results.basin1.z == d);
    %         O2_measured = res.T(res.depth1 == d);
    %         day_measured = res.date(res.depth1 == d);
    %         day_measured = day_measured(~isnan(O2_measured));
    %         O2_measured = O2_measured(~isnan(O2_measured));

    %         O2_mod = MyLake_results.basin1.concentrations.O2(zinx,:)'/1000;
    %         [T_date,loc_sim, loc_obs] = intersect(MyLake_results.basin1.days, day_measured);

    %         % rmsd_O2 = rmsd_O2 + RMSE(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
    %         rmsd_O2 = rmsd_O2 + sqrt(mean((O2_mod(loc_sim, 1)-O2_measured(loc_obs, 1)).^2));
    %     end

    %     zinx=find(MyLake_results.basin1.z<4);
    %     TP_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:)+MyLake_results.basin1.concentrations.PP(zinx,:) + MyLake_results.basin1.concentrations.DOP(zinx,:) + MyLake_results.basin1.concentrations.POP(zinx,:))', 2);
    %     Chl_mod = mean((MyLake_results.basin1.concentrations.Chl(zinx,:)+MyLake_results.basin1.concentrations.C(zinx,:))', 2);
    %     P_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:))', 2);
    %     POP_mod = mean((MyLake_results.basin1.concentrations.POP(zinx,:) + MyLake_results.basin1.concentrations.PP(zinx,:))', 2);

    %     load 'obs/store_obs/TOTP.dat' % measured
    %     load 'obs/store_obs/Cha.dat' % measured
    %     load 'obs/store_obs/PO4.dat' % measured
    %     load 'obs/store_obs/Part.dat' % measured


    %     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, TOTP(:,1)));
    %     rmsd_TOTP = sqrt(mean((TP_mod(loc_sim, 1)-TOTP(loc_obs, 2)).^2));


    %     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha(:,1)));
    %     rmsd_Chl = sqrt(mean((Chl_mod(loc_sim, 1)-Cha(loc_obs, 2)).^2));


    %     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    %     rmsd_PO4 = sqrt(mean((P_mod(loc_sim, 1)-PO4(loc_obs, 2)).^2));


    %     [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    %     rmsd_PP = sqrt(mean((POP_mod(loc_sim, 1)-Part(loc_obs, 2)).^2));


    %     disp('RMSD 3xRMSE(P)+RMSE(O2):')
    %     disp(sum([3*rmsd_TOTP, 3*rmsd_Chl, 3*rmsd_PO4, 3*rmsd_PP, rmsd_O2]))
    %     disp('RMSD = RMSE(P)+RMSE(O2):')
    %     disp(sum([rmsd_TOTP, rmsd_Chl, rmsd_PO4, rmsd_PP, rmsd_O2]))
    % end


    % toc



end
% end
