function x = calibration_routine_sed()
tic
format shortEng
format compact
% parpool
% gaoptions = optimoptions('ga','UseParallel',true);

% x = [0.032; 0.3627];
% lb = x0*0.1;
% ub = x0*10;

[lake_params, sediment_params] = load_params();


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
sediment_params{75} = 0.6;%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo
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
lake_params{19} = 150;              % 19    scaling factor for inflow concentration of POC (-)


lake_params{18} = 10; % 1.5525e+000   % 18    Isc C
lake_params{22} = 10; % 1.5525e+000   % 22    Isc Chl

sediment_params{8} = 300;%    Km FeOH3
sediment_params{9} = 300;%    Km FeOOH
% ======================================================================


% Initial values
% Niva sediment cores & inputs scaled & k_chl=3; err= r^2*RMSD, res=~850.34  % ====================================================
x(1) = lake_params{52}; % 65.1237e+000   %    accel
x(2) = sediment_params{1};  %   'k_Chl',                 %        % 1
x(3) = sediment_params{2};  %  'k_POP',                 %        % 1
x(4) = sediment_params{3};  % 'k_POC',                  %        % 0.01
x(5) = sediment_params{4};  %  'k_DOP',                 %        % 1
x(6) = sediment_params{5};  % 'k_DOC',                  %        % 1
x(7) = sediment_params{23};  %     'k_pdesorb_a',         %
x(8) = sediment_params{24};  %     'k_pdesorb_b',         %
x(9) = sediment_params{54};  %     'k_pdesorb_c',         %

% SO4 boundary
x(10) = sediment_params{75};%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo

% for cores too (scaling unknown inputs):
x(11) = lake_params{25};%    Scaling factor for inflow concentration of O2 (-)
x(12) = lake_params{27};%    Scaling factor for inflow concentration of NO3 (-)
x(13) = lake_params{34};%    Scaling factor for inflow concentration of Fe3 (-)
x(14) = lake_params{35};%    Scaling factor for inflow concentration of Al3 (-)
x(15) = lake_params{37};%    Scaling factor for inflow concentration of CaCO3 (-)

% P minerals:
x(16) = sediment_params{31};%    k_apa_pre
x(17) = sediment_params{32};%    k_apa_pre
x(18) = sediment_params{40};%    k_viv_pre
x(19) = sediment_params{41};%    k_viv_pre

x(20) = sediment_params{8};%    Km FeOH3
x(21) = sediment_params{9};%    Km FeOOH


lb = [1; 0.01; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001;   0;    0;    0;   0;   0;   0; 1e-6; 1e-6; 1e-6; 1e-6; 1; 1];
ub = [300;  1;   0.1;   0.1;   0.1;   0.1;  1000; 1000;   1000; 100;    2;  1000; 1000; 1000; 1000;   10;   10;   10;   10; 1000; 1000];


fcns = {@gaplotscorediversity, @gaplotstopping, @gaplotgenealogy, @gaplotscores, @gaplotdistance, @gaplotselection, @gaplotmaxconstr, @gaplotbestf, @gaplotbestindiv, @gaplotexpectation, @gaplotrange, @gaplotpareto, @gaplotparetodistance, @gaplotrankhist, @gaplotspread};

population_size = 72;  % Populations size for each generation of the genetic algorithm
max_generations = 7;  % How many generations to run the genetic algorithm for
parallelize     = true;

% options = gaoptimset('Display','iter','UseParallel', true, 'TolFun', 1e-2, 'PlotFcns', fcns);
options = optimoptions('ga', 'MaxGenerations', max_generations, 'PopulationSize', population_size, 'UseParallel', parallelize);
x = ga(@opt_fun,length(x),[],[],[],[],lb,ub, @nonlcon, options)




%% opt_fun: function which we are going to minimize
function [res] = opt_fun(x)

[lake_params, sediment_params] = load_params();


% new added for cores
lake_params{52} = x(1); % 65.1237e+000   %    accel
sediment_params{1} = x(2);  %   'k_Chl',                 %        % 1
sediment_params{2} = x(3);  %  'k_POP',                 %        % 1
sediment_params{3} = x(4);  % 'k_POC',                  %        % 0.01
sediment_params{4} = x(5);  %  'k_DOP',                 %        % 1
sediment_params{5} = x(6);  % 'k_DOC',                  %        % 1
sediment_params{23} = x(7);  %     'k_pdesorb_a',         %
sediment_params{24} = x(8);  %     'k_pdesorb_b',         %
sediment_params{54} = x(9);  %     'k_pdesorb_c',         %

% SO4 boundary
sediment_params{75} = x(10);%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo

% for cores too (scaling unknown inputs):
lake_params{25} = x(11);%    Scaling factor for inflow concentration of O2 (-)
lake_params{27} = x(12);%    Scaling factor for inflow concentration of NO3 (-)
lake_params{34} = x(13);%    Scaling factor for inflow concentration of Fe3 (-)
lake_params{35} = x(14);%    Scaling factor for inflow concentration of Al3 (-)
lake_params{37} = x(15);%    Scaling factor for inflow concentration of CaCO3 (-)

% P minerals:
sediment_params{31} = x(16);%    k_apa_pre
sediment_params{32} = x(17);%    k_apa_pre
sediment_params{40} = x(18);%    k_viv_pre
sediment_params{41} = x(19);%    k_viv_pre

sediment_params{8} = x(20);%    Km FeOH3
sediment_params{9} = x(21);%    Km FeOOH


% modifications:
sediment_params{73} = 24;
sediment_params{74} = 0; % pH module off, const pH = 8

run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID;
m_start=[2000, 1, 1]; % Do not change this date if you are calibrating the cores (using relative dates in the code) or check it
m_stop=[2013, 10, 31]; %
run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files
is_save_results = false;


% disp(datetime('now'));

try
    name_of_scenario = 'IO/store_INCAP_input_baseline_mod.txt';
    [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, name_of_scenario, use_INCA, run_INCA, run_ID, clim_ID, is_save_results); % runs the model and outputs obs and sim


    load('Postproc_code/Vansjo/VAN1_data_2017_02_28_10_55.mat')

    depths = [5;10;15;20;25;30;35;40];
    nrmsd_mean_O2 = 0;
    rsquared_O2 = 0;


    for i=1:size(depths,1)
        d = depths(i);
        zinx=find(MyLake_results.basin1.z == d);
        O2_measured = res.T(res.depth1 == d);
        day_measured = res.date(res.depth1 == d);
        day_measured = day_measured(~isnan(O2_measured));
        O2_measured = O2_measured(~isnan(O2_measured));

        O2_mod = MyLake_results.basin1.concentrations.O2(zinx,:)'/1000;
        [T_date,loc_sim, loc_obs] = intersect(MyLake_results.basin1.days, day_measured);

        nrmsd_mean_O2(i) = nrmsd_mean(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        rsquared_O2(i) = rsquared(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        % nrmsd_mean_O2 = nrmsd_mean_O2 + sqrt(mean((O2_mod(loc_sim, 1)-O2_measured(loc_obs, 1)).^2));
    end

    % P forms measured in water-column
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
    nrmsd_mean_TOTP = nrmsd_mean(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));
    rsquared_TOTP = rsquared(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha(:,1)));
    nrmsd_mean_Chl = nrmsd_mean(Chl_mod(loc_sim, 1), Cha(loc_obs, 2));
    rsquared_Chl = rsquared(Chl_mod(loc_sim, 1), Cha(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    nrmsd_mean_PO4 = nrmsd_mean(P_mod(loc_sim, 1), PO4(loc_obs, 2));
    rsquared_PO4 = rsquared(P_mod(loc_sim, 1), PO4(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    nrmsd_mean_PP = nrmsd_mean(POP_mod(loc_sim, 1), Part(loc_obs, 2));
    rsquared_PP = rsquared(POP_mod(loc_sim, 1), Part(loc_obs, 2));

    % Sediment cores measured in October 2013
    load 'obs/store_obs/P_HCl_sed.dat'
    load 'obs/store_obs/P_Ca_sed.dat' %
    load 'obs/store_obs/P_Org_sed.dat'
    load 'obs/store_obs/P_Al_sed.dat'
    load 'obs/store_obs/P_Fe_sed.dat'  %
    load 'obs/store_obs/P_H2O_sed.dat'
    load 'obs/store_obs/S_sed.dat' %
    load 'obs/store_obs/Fe2_sed.dat' %
    load 'obs/store_obs/Ca_sed.dat' %
    load 'obs/store_obs/Al3_sed.dat'
    load 'obs/store_obs/PO4_sed.dat' %


    sed_core_date = 735523; % =datenum('14-Oct-2013','dd-mmm-yyyy')

    [~,idx_date_sed_cores,~] = intersect(MyLake_results.basin1.days, 735523);
    idx_depthx_sed_cores =  floor(PO4_sed(:,1)/Sediment_results.basin1.z(end)*(Sediment_results.basin1.params.n-1));

    nrmsd_mean_PO4_sed = nrmsd_mean(30.973*Sediment_results.basin1.concentrations.PO4(idx_depthx_sed_cores,idx_date_sed_cores), PO4_sed(:,2)+P_H2O_sed(:,2));
    rsquared_PO4_sed = rsquared(30.973*Sediment_results.basin1.concentrations.PO4(idx_depthx_sed_cores,idx_date_sed_cores), PO4_sed(:,2)+P_H2O_sed(:,2));

    nrmsd_mean_Ca_sed = nrmsd_mean(40.0784*Sediment_results.basin1.concentrations.Ca2(idx_depthx_sed_cores,idx_date_sed_cores), Ca_sed(:,2));
    rsquared_Ca_sed = rsquared(40.0784*Sediment_results.basin1.concentrations.Ca2(idx_depthx_sed_cores,idx_date_sed_cores), Ca_sed(:,2));

    nrmsd_mean_Fe_sed = nrmsd_mean(55.8452*Sediment_results.basin1.concentrations.Fe2(idx_depthx_sed_cores,idx_date_sed_cores), Fe2_sed(:,2));
    rsquared_Fe_sed = rsquared(55.8452*Sediment_results.basin1.concentrations.Fe2(idx_depthx_sed_cores,idx_date_sed_cores), Fe2_sed(:,2));

    nrmsd_mean_S_sed = nrmsd_mean(32.0655*Sediment_results.basin1.concentrations.SO4(idx_depthx_sed_cores,idx_date_sed_cores), S_sed(:,2));
    rsquared_S_sed = rsquared(32.0655*Sediment_results.basin1.concentrations.SO4(idx_depthx_sed_cores,idx_date_sed_cores), S_sed(:,2));

    nrmsd_mean_P_Fe_sed = nrmsd_mean(...
        30.973*Sediment_results.basin1.concentrations.PO4adsa(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        30.973*Sediment_results.basin1.concentrations.PO4adsb(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        2*30.973*Sediment_results.basin1.concentrations.Fe3PO42(idx_depthx_sed_cores,idx_date_sed_cores), ...
        P_Fe_sed(:,2));
    rsquared_P_Fe_sed = rsquared(...
        30.973*Sediment_results.basin1.concentrations.PO4adsa(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        30.973*Sediment_results.basin1.concentrations.PO4adsb(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        2*30.973*Sediment_results.basin1.concentrations.Fe3PO42(idx_depthx_sed_cores,idx_date_sed_cores), ...
        P_Fe_sed(:,2));

    nrmsd_mean_P_Ca_sed = nrmsd_mean(2*30.973*Sediment_results.basin1.concentrations.Ca3PO42(idx_depthx_sed_cores,idx_date_sed_cores), P_Ca_sed(:,2));
    rsquared_P_Ca_sed = rsquared(2*30.973*Sediment_results.basin1.concentrations.Ca3PO42(idx_depthx_sed_cores,idx_date_sed_cores), P_Ca_sed(:,2));

    % nrmsd_mean_POP_sed = nrmsd_mean(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.POP(idx_depthx_sed_cores,idx_date_sed_cores), P_Org_sed(:,2));
    % rsquared_POP_sed = rsquared(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.POP(idx_depthx_sed_cores,idx_date_sed_cores), P_Org_sed(:,2));

    nrmsd_mean_P_Al_sed = nrmsd_mean(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.PO4adsc(idx_depthx_sed_cores,idx_date_sed_cores), P_Al_sed(:,2));
    rsquared_P_Al_sed = nrmsd_mean(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.PO4adsc(idx_depthx_sed_cores,idx_date_sed_cores), P_Al_sed(:,2));

    x'

    % res = sum([3*nrmsd_mean_TOTP, 3*nrmsd_mean_Chl, 3*nrmsd_mean_PO4, 3*nrmsd_mean_PP, nrmsd_mean_O2])
    % res = sum([- (rsquared_TOTP - 1), - (rsquared_Chl - 1), - (rsquared_PO4 - 1), - (rsquared_PP - 1), mean(- (rsquared_O2 + 1))])

    % just nrmsd_mean
    res = sum([3*nrmsd_mean_PO4_sed, nrmsd_mean_Ca_sed, 3*nrmsd_mean_Fe_sed, 3*nrmsd_mean_S_sed, 3*nrmsd_mean_P_Fe_sed, nrmsd_mean_P_Ca_sed, nrmsd_mean_P_Al_sed])





catch ME
    fprintf('\tID: %s\n', ME.identifier)
    fprintf('\tMessage: %s\n', ME.message)
    fprintf('\tStack::\n')
    for k=1:length(ME.stack)
        disp(ME.stack(k))
    end
    res = NaN
end



function [c,ceq] = nonlcon(x)
c = -x;
ceq = [];


