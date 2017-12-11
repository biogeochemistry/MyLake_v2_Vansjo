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

% Niva sediment cores & inputs scaled & k_chl=3; err= r^2*RMSD, res=~850.34  % ====================================================
x = [50; 0.111452954004814; 0.0756702079558015; 0.0633637582049403; 0.100000000000000; 0.0957726460081909; 1.96425357122318; 19.1927004493911; 19.1927004493911; 0.1; 1; 1; 1; 1; 1; 1];



lb = [1; 0.01; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001;   0;    0; 0;    0;   0;   0;   0];
ub = [100;  1;   0.1;   0.1;   0.1;   0.1;   100;   100;   100; 100;  100; 2;  100; 100; 100; 100];


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



run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID;
m_start=[2000, 1, 1]; % Do not change this date if you are calibrating the cores (using relative dates in the code)
m_stop=[2013, 10, 31]; %
run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files
is_save_results = false;


% disp(datetime('now'));

try

    [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, is_save_results); % runs the model and outputs obs and sim


    load('Postproc_code/Vansjo/VAN1_data_2017_02_28_10_55.mat')

    depths = [5;10;15;20;25;30;35;40];
    rmsd_O2 = 0;
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

        rmsd_O2(i) = rmsd(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        rsquared_O2(i) = rsquared(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        % rmsd_O2 = rmsd_O2 + sqrt(mean((O2_mod(loc_sim, 1)-O2_measured(loc_obs, 1)).^2));
    end

    % P forms measured in water-column
    zinx=find(MyLake_results.basin1.z<4);
    TP_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:)+MyLake_results.basin1.concentrations.PP(zinx,:) + MyLake_results.basin1.concentrations.DOP(zinx,:) + MyLake_results.basin1.concentrations.POP(zinx,:))', 2);
    Chl_mod = mean((MyLake_results.basin1.concentrations.Chl(zinx,:)+MyLake_results.basin1.concentrations.C(zinx,:))', 2);
    P_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:))', 2);
    POP_mod = mean((MyLake_results.basin1.concentrations.POP(zinx,:) + MyLake_results.basin1.concentrations.PP(zinx,:))', 2);

    load 'obs/store_obs/TOTP.dat' % measured
    % load 'obs/store_obs/Cha.dat' % measured
    load 'obs/store_obs/Cha_aquaM_march_2017.dat' % measured
    load 'obs/store_obs/PO4.dat' % measured
    load 'obs/store_obs/Part.dat' % measured


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, TOTP(:,1)));
    rmsd_TOTP = rmsd(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));
    rsquared_TOTP = rsquared(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha_aquaM_march_2017(:,1)));
    rmsd_Chl = rmsd(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));
    rsquared_Chl = rsquared(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    rmsd_PO4 = rmsd(P_mod(loc_sim, 1), PO4(loc_obs, 2));
    rsquared_PO4 = rsquared(P_mod(loc_sim, 1), PO4(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    rmsd_PP = rmsd(POP_mod(loc_sim, 1), Part(loc_obs, 2));
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

    rmsd_PO4_sed = rmsd(30.973*Sediment_results.basin1.concentrations.PO4(idx_depthx_sed_cores,idx_date_sed_cores), PO4_sed(:,2)+P_H2O_sed(:,2));
    rsquared_PO4_sed = rsquared(30.973*Sediment_results.basin1.concentrations.PO4(idx_depthx_sed_cores,idx_date_sed_cores), PO4_sed(:,2)+P_H2O_sed(:,2));

    rmsd_Ca_sed = rmsd(40.0784*Sediment_results.basin1.concentrations.Ca2(idx_depthx_sed_cores,idx_date_sed_cores), Ca_sed(:,2));
    rsquared_Ca_sed = rsquared(40.0784*Sediment_results.basin1.concentrations.Ca2(idx_depthx_sed_cores,idx_date_sed_cores), Ca_sed(:,2));

    rmsd_Fe_sed = rmsd(55.8452*Sediment_results.basin1.concentrations.Fe2(idx_depthx_sed_cores,idx_date_sed_cores), Fe2_sed(:,2));
    rsquared_Fe_sed = rsquared(55.8452*Sediment_results.basin1.concentrations.Fe2(idx_depthx_sed_cores,idx_date_sed_cores), Fe2_sed(:,2));

    rmsd_S_sed = rmsd(32.0655*Sediment_results.basin1.concentrations.SO4(idx_depthx_sed_cores,idx_date_sed_cores), S_sed(:,2));
    rsquared_S_sed = rsquared(32.0655*Sediment_results.basin1.concentrations.SO4(idx_depthx_sed_cores,idx_date_sed_cores), S_sed(:,2));

    rmsd_P_Fe_sed = rmsd(...
        30.973*Sediment_results.basin1.concentrations.PO4adsa(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        30.973*Sediment_results.basin1.concentrations.PO4adsb(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        2*30.973*Sediment_results.basin1.concentrations.Fe3PO42(idx_depthx_sed_cores,idx_date_sed_cores), ...
        P_Fe_sed(:,2));
    rsquared_P_Fe_sed = rsquared(...
        30.973*Sediment_results.basin1.concentrations.PO4adsa(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        30.973*Sediment_results.basin1.concentrations.PO4adsb(idx_depthx_sed_cores,idx_date_sed_cores) + ...
        2*30.973*Sediment_results.basin1.concentrations.Fe3PO42(idx_depthx_sed_cores,idx_date_sed_cores), ...
        P_Fe_sed(:,2));

    rmsd_P_Ca_sed = rmsd(2*30.973*Sediment_results.basin1.concentrations.Ca3PO42(idx_depthx_sed_cores,idx_date_sed_cores), P_Ca_sed(:,2));
    rsquared_P_Ca_sed = rsquared(2*30.973*Sediment_results.basin1.concentrations.Ca3PO42(idx_depthx_sed_cores,idx_date_sed_cores), P_Ca_sed(:,2));

    rmsd_POP_sed = rmsd(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.POP(idx_depthx_sed_cores,idx_date_sed_cores), P_Org_sed(:,2));
    rsquared_POP_sed = rsquared(30.973*Sediment_results.basin1.params.Pz1*Sediment_results.basin1.concentrations.POP(idx_depthx_sed_cores,idx_date_sed_cores), P_Org_sed(:,2));

    x'

    % res = sum([3*rmsd_TOTP, 3*rmsd_Chl, 3*rmsd_PO4, 3*rmsd_PP, rmsd_O2])
    % res = sum([- (rsquared_TOTP - 1), - (rsquared_Chl - 1), - (rsquared_PO4 - 1), - (rsquared_PP - 1), mean(- (rsquared_O2 + 1))])

    k_chl = 3;

    res = sum([- (rsquared_TOTP - 1) .* rmsd_TOTP, - (rsquared_Chl - 1) .* rmsd_Chl * k_chl, - (rsquared_PO4 - 1) .* rmsd_PO4, - (rsquared_PP - 1) .* rmsd_PP, mean(- (rsquared_O2 - 1) .* rmsd_O2), - (rsquared_PO4_sed - 1) .* rmsd_PO4_sed, - (rsquared_Ca_sed - 1) .* rmsd_Ca_sed, - (rsquared_Fe_sed - 1) .* rmsd_Fe_sed, - (rsquared_S_sed - 1) .* rmsd_S_sed, - (rsquared_P_Fe_sed - 1) .* rmsd_P_Fe_sed, - (rsquared_P_Ca_sed - 1) .* rmsd_P_Ca_sed, - (rsquared_POP_sed - 1) .* rmsd_POP_sed])





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


