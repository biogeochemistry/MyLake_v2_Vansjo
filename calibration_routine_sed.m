function x = calibration_routine()
tic
format shortEng
format compact
% parpool
% gaoptions = optimoptions('ga','UseParallel',true);

% x = [0.032; 0.3627];
% lb = x0*0.1;
% ub = x0*10;

[lake_params, sediment_params] = load_params();

% NIVA latest RMSD only Chl result (++++Best so far++++):
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


% Niva sediment cores & inputs scaled & k_chl=3; err= r^2*RMSD, res=~850.34  % ====================================================
x = [50; 0.111452954004814; 0.0756702079558015; 0.0633637582049403; 0.100000000000000; 0.0957726460081909; 1.96425357122318; 19.1927004493911; 19.1927004493911; 0.1; 1; 1; 1; 1; 1; 1; 0.00037; 0.37; 0.00037; 0.37; 100/2.5; 100/2.5];



lb = [1; 0.01; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001;   0;    0; 0;    0;   0;   0;   0; 1e-6; 1e-6; 1e-6; 1e-6; 1; 1];
ub = [100;  1;   0.1;   0.1;   0.1;   0.1;   100;   100;   100; 100;  100; 2;  100; 100; 100; 100;   10;   10;   10;   10; 1000; 1000];



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
sediment_params{1} = x(1); % 65.1237e+000   %    accel
sediment_params{1} = x(2);  %   'k_Chl',                 %        % 1
sediment_params{2} = x(3);  %  'k_POP',                 %        % 1
sediment_params{3} = x(4);  % 'k_POC',                  %        % 0.01
sediment_params{4} = x(5);  %  'k_DOP',                 %        % 1
sediment_params{5} = x(6);  % 'k_DOC',                  %        % 1
sediment_params{23} = x(7);  %     'k_pdesorb_a',         %
sediment_params{24} = x(8);  %     'k_pdesorb_b',         %
sediment_params{54} = x(9);  %     'k_pdesorb_c',         %

% SO4 boundary
lake_params{75} = x(10);%    % flux of sulphate from bottom of the sediment. Custom boundary condition for Vansjo

% for cores too (scaling unknown inputs):
lake_params{22} = x(11);%    scaling factor for inflow concentration of Chl a (-)
lake_params{25} = x(12);%    Scaling factor for inflow concentration of O2 (-)
lake_params{27} = x(13);%    Scaling factor for inflow concentration of NO3 (-)
lake_params{34} = x(14);%    Scaling factor for inflow concentration of Fe3 (-)
lake_params{35} = x(15);%    Scaling factor for inflow concentration of Al3 (-)
lake_params{37} = x(16);%    Scaling factor for inflow concentration of CaCO3 (-)

% P minerals:
lake_params{31} = x(17);%    k_apa_pre
lake_params{32} = x(18);%    k_apa_pre
lake_params{40} = x(19);%    k_viv_pre
lake_params{41} = x(20);%    k_viv_pre

sediment_params{8} = x(21);%    Km FeOH3
sediment_params{9} = x(22);%    Km FeOOH


% modifications:
sediment_params{73} = 48;
sediment_params{74} = 0; % pH module off, const pH = 8

run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID;
m_start=[2000, 1, 1]; % Do not change this date if you are calibrating the cores (using relative dates in the code)
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
    nrmsd_O2 = 0;
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

        nrmsd_O2(i) = nrmsd(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        rsquared_O2(i) = rsquared(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        % nrmsd_O2 = nrmsd_O2 + sqrt(mean((O2_mod(loc_sim, 1)-O2_measured(loc_obs, 1)).^2));
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
    nrmsd_TOTP = nrmsd(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));
    rsquared_TOTP = rsquared(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha(:,1)));
    nrmsd_Chl = nrmsd(Chl_mod(loc_sim, 1), Cha(loc_obs, 2));
    rsquared_Chl = rsquared(Chl_mod(loc_sim, 1), Cha(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    nrmsd_PO4 = nrmsd(P_mod(loc_sim, 1), PO4(loc_obs, 2));
    rsquared_PO4 = rsquared(P_mod(loc_sim, 1), PO4(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    nrmsd_PP = nrmsd(POP_mod(loc_sim, 1), Part(loc_obs, 2));
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

    nrmsd_PO4_sed = nrmsd(30.973*Sediment_results.basin1.concentrations.PO4(idx_depthx_sed_cores,idx_date_sed_cores), PO4_sed(:,2)+P_H2O_sed(:,2));
    rsquared_PO4_sed = rsquared(30.973*Sediment_results.basin1.concentrations.PO4(idx_depthx_sed_cores,idx_date_sed_cores), PO4_sed(:,2)+P_H2O_sed(:,2));

    nrmsd_Ca_sed = nrmsd(40.0784*Sediment_results.basin1.concentrations.Ca2(idx_depthx_sed_cores,idx_date_sed_cores), Ca_sed(:,2));
    rsquared_Ca_sed = rsquared(40.0784*Sediment_results.basin1.concentrations.Ca2(idx_depthx_sed_cores,idx_date_sed_cores), Ca_sed(:,2));

    nrmsd_Fe_sed = nrmsd(55.8452*Sediment_results.basin1.concentrations.Fe2(idx_depthx_sed_cores,idx_date_sed_cores), Fe2_sed(:,2));
    rsquared_Fe_sed = rsquared(55.8452*Sediment_results.basin1.concentrations.Fe2(idx_depthx_sed_cores,idx_date_sed_cores), Fe2_sed(:,2));

    nrmsd_S_sed = nrmsd(32.0655*Sediment_results.basin1.concentrations.SO4(idx_depthx_sed_cores,idx_date_sed_cores), S_sed(:,2));
    rsquared_S_sed = rsquared(32.0655*Sediment_results.basin1.concentrations.SO4(idx_depthx_sed_cores,idx_date_sed_cores), S_sed(:,2));

    nrmsd_P_Fe_sed = nrmsd(...
        30.973*Sediment_results.basin1.concentrations.PO4adsa(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        30.973*Sediment_results.basin1.concentrations.PO4adsb(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        2*30.973*Sediment_results.basin1.concentrations.Fe3PO42(idx_depthx_sed_cores,idx_date_sed_cores), ...
        P_Fe_sed(:,2));
    rsquared_P_Fe_sed = rsquared(...
        30.973*Sediment_results.basin1.concentrations.PO4adsa(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        30.973*Sediment_results.basin1.concentrations.PO4adsb(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        2*30.973*Sediment_results.basin1.concentrations.Fe3PO42(idx_depthx_sed_cores,idx_date_sed_cores), ...
        P_Fe_sed(:,2));

    nrmsd_P_Ca_sed = nrmsd(2*30.973*Sediment_results.basin1.concentrations.Ca3PO42(idx_depthx_sed_cores,idx_date_sed_cores), P_Ca_sed(:,2));
    rsquared_P_Ca_sed = rsquared(2*30.973*Sediment_results.basin1.concentrations.Ca3PO42(idx_depthx_sed_cores,idx_date_sed_cores), P_Ca_sed(:,2));

    % nrmsd_POP_sed = nrmsd(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.POP(idx_depthx_sed_cores,idx_date_sed_cores), P_Org_sed(:,2));
    % rsquared_POP_sed = rsquared(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.POP(idx_depthx_sed_cores,idx_date_sed_cores), P_Org_sed(:,2));

    nrmsd_P_Al_sed = nrmsd(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.PO4adsc(idx_depthx_sed_cores,idx_date_sed_cores), P_Al_sed(:,2));
    rsquared_P_Al_sed = nrmsd(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.PO4adsc(idx_depthx_sed_cores,idx_date_sed_cores), P_Al_sed(:,2));

    x'

    % res = sum([3*nrmsd_TOTP, 3*nrmsd_Chl, 3*nrmsd_PO4, 3*nrmsd_PP, nrmsd_O2])
    % res = sum([- (rsquared_TOTP - 1), - (rsquared_Chl - 1), - (rsquared_PO4 - 1), - (rsquared_PP - 1), mean(- (rsquared_O2 + 1))])

    % just nrmsd
    res = sum([nrmsd_TOTP, nrmsd_Chl, nrmsd_PO4, nrmsd_PP, mean(nrmsd_O2), 2*nrmsd_PO4_sed, 2*nrmsd_Ca_sed, 2*nrmsd_Fe_sed, nrmsd_S_sed, 2*nrmsd_P_Fe_sed, 2*nrmsd_P_Ca_sed, nrmsd_P_Al_sed])





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


