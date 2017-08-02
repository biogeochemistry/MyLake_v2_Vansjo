function x = calibration_routine()
tic
format shortEng
format compact
% parpool
% gaoptions = optimoptions('ga','UseParallel',true);

% x = [0.032; 0.3627];
% lb = x0*0.1;
% ub = x0*10;

x = [472.4167e-003,   192.8083e-003,     1.0000e+000,   806.4479e-003,    59.1485e-003,   100.0000e-003,     1.0000e+000,     1.8765e+000,    10.0000e-003,     2.0000e+000,     2.1591e+000,     0.0000e+000,     0.0000e+000,   148.5501e-003,   100.0000e-006,    10.0000e-006,     5.0771e-003,    38.0006e-003]; % RMSD 135.6022

x(12) = 1; % 23    scaling factor for inflow concentration of POC  (-)
x(9) = 0.01;  % % 8  settling velocity for S (m day-1)
x(5) = 0.01; % 18    Settling velocity for Chl2 a (m day-1)

lb = [0.05, 0.1, 1, 0.2, 0.05, 0.1, 1, 0.2, 0.01, 1, 1,  0, 0, 0, 1e-5, 1e-5, 0.005, 0.005];
ub = [0.5, 0.3, 1.5, 2, 0.5, 0.3, 1.5, 2, 1, 1e5, 100,  10, 100, 100, 1e-4, 1e-4, 0.045, 0.045];


fcns = {@gaplotscorediversity, @gaplotstopping, @gaplotgenealogy, @gaplotscores, @gaplotdistance, @gaplotselection, @gaplotmaxconstr, @gaplotbestf, @gaplotbestindiv, @gaplotexpectation, @gaplotrange, @gaplotpareto, @gaplotparetodistance, @gaplotrankhist, @gaplotspread};

population_size = 72;  % Populations size for each generation of the genetic algorithm
max_generations = 7;  % How many generations to run the genetic algorithm for
parallelize     = true; % 15 generation takes 12 hours on 24 cores

% options = gaoptimset('Display','iter','UseParallel', true, 'TolFun', 1e-2, 'PlotFcns', fcns);
options = optimoptions('ga', 'MaxGenerations', max_generations, 'PopulationSize', population_size, 'UseParallel', parallelize);

x = ga(@opt_fun,length(x),[],[],[],[],lb,ub, @nonlcon, options)




%% opt_fun: function which we are going to minimize
function [res] = opt_fun(x)

[lake_params, sediment_params] = load_params();




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


lake_params{28 -5} = x(12); % 23    scaling factor for inflow concentration of DOC  (-)
lake_params{24 -5} = x(13); % 19    scaling factor for inflow concentration of POC (-)
lake_params{39 -5} = x(14); % 34    Scaling factor for inflow concentration of Fe3 (-)

lake_params{15 -5} = x(15); % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{59 -5} = x(16); % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{17 -5} = x(17); % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{60 -5} = x(18); % 17    Optical cross_section of chlorophyll (m2 mg-1)




run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID
m_start=[2002, 1, 1]; %
m_stop=[2008, 12, 31]; %
run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files
is_save_results = false;

% sediment_params{56} = 1; % Only for preliminary calibration: coarse time step for chemical and sediment modules


disp(datetime('now'));

try

    [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, is_save_results); % runs the model and outputs obs and sim


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

        rmsd_O2 = rmsd_O2 + RMSE(O2_mod(loc_sim, 1), O2_measured(loc_obs, 1));
        % rmsd_O2 = rmsd_O2 + sqrt(mean((O2_mod(loc_sim, 1)-O2_measured(loc_obs, 1)).^2));
    end

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
    rmsd_TOTP = RMSE(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha_aquaM_march_2017(:,1)));
    rmsd_Chl = RMSE(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    rmsd_PO4 = RMSE(P_mod(loc_sim, 1), PO4(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    rmsd_PP = RMSE(POP_mod(loc_sim, 1), Part(loc_obs, 2));


    x'

    res = sum([3*rmsd_TOTP, 5*rmsd_Chl, 5*rmsd_PO4, 3*rmsd_PP, rmsd_O2])

catch ME
    fprintf('\tID: %s\n', ME.identifier)
    fprintf('\tMessage: %s\n', ME.message)
    fprintf('\tStack::\n')
    for k=1:length(ME.stack)
        disp(ME.stack(k))
    end
    res = NaN
end


function r = RMSE(y, yhat)
    r = sqrt(mean((y-yhat).^2));

function [c,ceq] = nonlcon(x)
c = -x;
ceq = [];
