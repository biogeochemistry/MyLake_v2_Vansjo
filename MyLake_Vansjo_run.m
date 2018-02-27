% for i = 1:1000

tic
disp('Started at:')
disp(datetime('now'));

is_metrics = true; % print metrics in the end

m_start=[2000, 1, 1]; %
m_stop=[2013, 12, 31]; %
% big_results = cell(1,no_runs);  % collects the results
% big_inputs = cell(1,no_runs);   % collects the inputs
save_initial_conditions = false; % save final concentrations as initial for the next run


[lake_params, sediment_params] = load_params();

name_of_scenario = 'IO/store_INCAP_input_baseline_mod.txt'
% name_of_scenario = 'IO/Scenarios/P_2016_cutoff.txt'
% file_name = 'IO/P_2016_cutoff_test.mat'
% name_of_scenario = 'IO/airT_Scenarios/T_only_RCP4_IPSL.txt'
% file_name = 'IO/airT_Scenarios/T_only_RCP4_IPSL.mat'


% inaccurate but faster:
sediment_params{73} = 96; %ts
sediment_params{74} = 0; %pH

% Niva res. NIVA: calibration of sed. with ts=24 and custom weights (3x,1x), nrmsd, pH 8, no POP, Km_Fe
% res ~= 21
file_name = 'IO/equal_rates_96.mat'


% new added for cores
sediment_params{1} = 1.0549e-01;  %   'k_Chl',                 %        % 1
sediment_params{2} = 1.2624e-02;  %  'k_POP',                 %        % 1
sediment_params{3} = 5.2341e-02;  % 'k_POC',                  %        % 0.01
sediment_params{4} = 1.2941e-02;  %  'k_DOP',                 %        % 1
sediment_params{5} = 8.7662e-02;  % 'k_DOC',                  %        % 1
sediment_params{23} = 6.3601e+00;  %     'k_pdesorb_a',         %
sediment_params{24} = 1.1171e+01;  %     'k_pdesorb_b',         %
sediment_params{54} = 4.9036e+01;  %     'k_pdesorb_c',         %

% SO4 boundary
sediment_params{75} = 9.1213e+01;%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo

% for cores too (scaling unknown inputs):
lake_params{22} = 6.9253e+01;%    scaling factor for inflow concentration of Chl a (-)
lake_params{25} = 8.4622e-01;%    Scaling factor for inflow concentration of O2 (-)
lake_params{27} = 1.6113e+01;%    Scaling factor for inflow concentration of NO3 (-)
lake_params{34} = 3.6356e+01;%    Scaling factor for inflow concentration of Fe3 (-)
lake_params{35} = 4.1063e+01;%    Scaling factor for inflow concentration of Al3 (-)
lake_params{37} = 6.4648e+01;%    Scaling factor for inflow concentration of CaCO3 (-)

% P minerals:
sediment_params{31} = 9.4278e-01;%    k_apa_pre
sediment_params{32} = 7.7780e+00;%    k_apa_pre
sediment_params{40} = 1.3434e+00;%    k_viv_pre
sediment_params{41} = 2.1799e+00;%    k_viv_pre

sediment_params{8} = 8.7728e+01;%    Km FeOH3
sediment_params{9} = 3.1972e+00;%    Km FeOOH

lake_params{24} = 1.0120e+00; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
lake_params{20} = 7.6529e-01; % 1   % 20    scaling factor for inflow concentration of TP (-)

lake_params{47} = 1.6558e-01; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{49} = 1.7861e-01; % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
lake_params{50} = 1.3772e+00; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{53} = 2.9236e-01; % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
lake_params{56} = 1.1681e-01; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} = 2.3063e-01; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
lake_params{58} = 1.4571e+00; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} = 3.2470e-01; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
lake_params{46} = 7.7897e-02; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{10} = 8.3890e-05; % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} = 7.3357e-05; % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{12} = 4.5000e-02; % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{55} = 4.5000e-02; % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% ======================================================================

% NOTE: trials
lake_params{24} = 1.0; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
lake_params{20} = 1.0;  % 20    scaling factor for inflow concentration of TP (-)
sediment_params{52} = 40;%    accel
lake_params{37} = 1;%    Scaling factor for inflow concentration of CaCO3 (-)
lake_params{50} = 2.5; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{58} = 2.5; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{47} = 0.05; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{56} = 0.05; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{46} = 0.05; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{53} = 10; % 638.9222e-003  % 53    Half saturation growth P level Chl1 (mg/m3)
lake_params{59} = 10; % 1.5525e+000   % 59    Half saturation growth P level Chl2 (mg/m3)

lake_params{10} = 1.03890e-05; % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} = 1.03357e-05; % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)


lake_params{18} = 10; % 1.5525e+000   % 18    Isc C
lake_params{22} = 10; % 1.5525e+000   % 22    Isc Chl


% Sediment cores:
lake_params{24} = 1.4; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
lake_params{20} = 1.9;  % 20    scaling factor for inflow concentration of TP (-)
lake_params{46} = 0.1; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{19} = 100+20;              % 19    scaling factor for inflow concentration of POC (-)


sediment_params{8} = 1.5*10*100*2.5;%    Km FeOH3 pvc 100, umol/g, rho = 2.5
sediment_params{9} = 1.5*10*100*2.5;%    Km FeOOH
sediment_params{15} = sediment_params{8};%   'Kin_FeOH3',         % 15       % the same as Km rho=2.5
sediment_params{16} = sediment_params{9};%   'Kin_FeOOH',         % 16       % the same as Km rho=2.5
lake_params{34} = 550; %    Scaling factor for inflow concentration of Fe3 (-)
lake_params{31} =  1; %  'I_scCa2',             % 31    Scaling factor for inflow concentration of Ca2 (-)
lake_params{37} = 5; % Isc CaCO3
% sediment_params{31} = 0.00037; %  'k_apa_pre',          % 31
% sediment_params{62} = 0; % 7.2; %   'alfa0',                % 62


sediment_params{23} = 40;  %     'k_pdesorb_a',         %
sediment_params{24} = 40;  %     'k_pdesorb_b',         %
sediment_params{54} = 40;  %     'k_pdesorb_c',         %
lake_params{35} = 0.001;%    Scaling factor for inflow concentration of Al3 (-)

% -> FeS -> FeS2 -> FeOOH
sediment_params{30} = 0.04;  %     'k_fe_pre',         %
sediment_params{45} = 0.12;  %   'k_FeSpre',         %
sediment_params{75} = 0.5; % 9.0;%    % flux of SO4  Vansjo
sediment_params{10} = 1000;%    Km SO4

% Apatite:
sediment_params{33} = 10^-10.22;
sediment_params{31} = 0.000037/3; %/10; % apa_pre
sediment_params{32} = 0.037; %  apa_dis


% % Vivenite
sediment_params{40} = 0.00037*10*1.5; % 0.00037; %,  'k_viv_pre',          % 40
sediment_params{41} = 0.37/7; % 0.37; %,  'k_viv_dis',             % 41
% sediment_params{40} = 0.00037*10*4; % 0.00037; %,  'k_viv_pre',          % 40
% sediment_params{41} = 0.2; % 0.37; %,  'k_viv_dis',             % 41

% FeCO3 and CaCO3
% sediment_params{37} = 180; %  'k_FeCO3_pre',        % 37      % Cappellen (1996)
% sediment_params{38} = 0.25; %     'k_FeCO3_dis',        % 38      % Cappellen (1996)
% sediment_params{34} = 0.04; %0.04,  'k_CaCO3_pre',        % 34      % Katsev (2013)
% sediment_params{35} = 0.05; %0.05,  'k_CaCO3_dis',           % 35      % Katsev (2013)



% NOTE: Chl changed here
lake_params{50} = 2; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{58} = 2; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{47} = 0.05; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{56} = 0.05; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{53} = 15; % 638.9222e-003  % 53    Half saturation growth P level Chl1 (mg/m3)
lake_params{59} = 15; % 1.5525e+000   % 59    Half saturation growth P level Chl2 (mg/m3)


% changed OM rates:
sediment_params{52} = 300;%    accel
sediment_params{1} = 1/13 * 0.1;  %   'k_Chl',                 %        % 1
sediment_params{2} = 1/13 * 0.012;  %  'k_POP',                 %        % 1
sediment_params{3} = 1/13 * 0.05;   % 'k_POC',                  %        % 0.01
sediment_params{4} = 1/13 * 0.013;  %  'k_DOP',                 %        % 1
sediment_params{5} = 1/13 * 0.088;  % 'k_DOC',                  %        % 1
% To play with this parameter
sediment_params{72} = 30; %     'effective_depth',     % 72           % depth below which the lake is affected by sediments, [m], if -1 (experimental) , then sediments below pycnocline


% Bioirrigation
sediment_params{62} = 7.2; %7.2,   'alfa0',                % 62


sediment_params{62} = 1; %0.15*0.5,   'Kd_fe2',           % 53



% rate of dissollution and precipitation should be equal for Viv and Apa
% Apatite:
sediment_params{33} = 10^-10.22;
sediment_params{31} = 0.000037/3/10*1.1; %/10; % apa_pre
sediment_params{32} = 0.037; %  apa_dis
% % Vivenite
sediment_params{40} = 0.00037*10*1.5*2; % 0.00037; %,  'k_viv_pre',          % 40
sediment_params{41} = 0.37/7/10; % 0.37; %,  'k_viv_dis',             % 41

lake_params{37} = 5/2; % Isc CaCO3
lake_params{34} = 600; %    Scaling factor for inflow concentration of Fe3 (-)
% To improove fit of Ca2=>PO4, we can play with initial CaCO3 concentration;

sediment_params{23} = 5;  %     'k_pdesorb_a',         %
sediment_params{24} = 5;  %     'k_pdesorb_b',         %
sediment_params{54} = 5;  %     'k_pdesorb_c',         %


% try
run_ID = 0;
clim_ID = 0;
run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, name_of_scenario, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim





disp('Saving results...')
save(file_name, 'MyLake_results', 'Sediment_results')
disp('Finished at:')
disp(datetime('now'));

if is_metrics == true

    load('Postproc_code/Vansjo/VAN1_data_2017_02_28_10_55.mat')

    depths = [5;10;15;20;25;30;35;40];
    rmsd_O2 = 0;


    for i=1:size(depths,1)
        d = depths(i);
        zinx=find(MyLake_results.basin1.z == d);
        O2_measured = res.T(res.depth1 == d);
        day_measured = res.date(res.depth1 == d);
        day_measured = day_measured(~isnan(O2_measured));
        O2_measured = O2_measured(~isnan(O2_measured));

        O2_mod = MyLake_results.basin1.concentrations.O2(zinx,:)'/1000;
        [T_date,loc_sim, loc_obs] = intersect(MyLake_results.basin1.days, day_measured);

        % rmsd_O2 = rmsd_O2 + RMSE(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        rmsd_O2 = rmsd_O2 + sqrt(mean((O2_mod(loc_sim, 1)-O2_measured(loc_obs, 1)).^2));
    end

    zinx=find(MyLake_results.basin1.z<4);
    TP_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:)+MyLake_results.basin1.concentrations.PP(zinx,:) + MyLake_results.basin1.concentrations.DOP(zinx,:) + MyLake_results.basin1.concentrations.POP(zinx,:))', 2);
    Chl_mod = mean((MyLake_results.basin1.concentrations.Chl(zinx,:)+MyLake_results.basin1.concentrations.C(zinx,:))', 2);
    P_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:))', 2);
    POP_mod = mean((MyLake_results.basin1.concentrations.POP(zinx,:) + MyLake_results.basin1.concentrations.PP(zinx,:))', 2);

    load 'obs/store_obs/TOTP.dat' % measured
    load 'obs/store_obs/Cha.dat' % measured
    load 'obs/store_obs/PO4.dat' % measured
    load 'obs/store_obs/Part.dat' % measured


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, TOTP(:,1)));
    rmsd_TOTP = sqrt(mean((TP_mod(loc_sim, 1)-TOTP(loc_obs, 2)).^2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha(:,1)));
    rmsd_Chl = sqrt(mean((Chl_mod(loc_sim, 1)-Cha(loc_obs, 2)).^2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    rmsd_PO4 = sqrt(mean((P_mod(loc_sim, 1)-PO4(loc_obs, 2)).^2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    rmsd_PP = sqrt(mean((POP_mod(loc_sim, 1)-Part(loc_obs, 2)).^2));


    disp('RMSD 3xRMSE(P)+RMSE(O2):')
    disp(sum([3*rmsd_TOTP, 3*rmsd_Chl, 3*rmsd_PO4, 3*rmsd_PP, rmsd_O2]))
    disp('RMSD = RMSE(P)+RMSE(O2):')
    disp(sum([rmsd_TOTP, rmsd_Chl, rmsd_PO4, rmsd_PP, rmsd_O2]))
end


toc




% end