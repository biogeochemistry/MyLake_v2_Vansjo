function x = calibration_routine()
tic
format shortEng
format compact
% parpool
% gaoptions = optimoptions('ga','UseParallel',true);
x = [0.1, 0.1368, 1.4647, 0.5293, 0.1, 0.1000, 1.3339, 0.8793, 0.05, 101.1213, 10];


% lb = x0*0.1;
% ub = x0*10;

lb = [0.05; 0.1; 1; 0.2; 0.05; 0.1; 1; 0.2; 0.01; 1; 1];
ub = [0.5; 0.3; 1.5; 2; 0.5; 0.3; 1.5; 2; 1; 1e6; 100];


fcns = {@gaplotscorediversity, @gaplotstopping, @gaplotgenealogy, @gaplotscores, @gaplotdistance, @gaplotselection, @gaplotmaxconstr, @gaplotbestf, @gaplotbestindiv, @gaplotexpectation, @gaplotrange, @gaplotpareto, @gaplotparetodistance, @gaplotrankhist, @gaplotspread};

population_size = 72;  % Populations size for each generation of the genetic algorithm
max_generations = 7;  % How many generations to run the genetic algorithm for
parallelize     = true; % 15 generation takes 12 hours on 24 cores

% options = gaoptimset('Display','iter','UseParallel', true, 'TolFun', 1e-2, 'PlotFcns', fcns);
options = optimoptions('ga', 'MaxGenerations', max_generations, 'PopulationSize', population_size, 'UseParallel', parallelize);

x = ga(@opt_fun,11,[],[],[],[],lb,ub, @nonlcon, options)

%% opt_fun: function which we are going to minimize
function [res] = opt_fun(x)

[lake_params, sediment_params] = load_params();


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files
is_save_results = false;

lake_params{52 -5} = x(1); % 9     settling velocity for Chl1 a (m day-1)
lake_params{54 -5} = x(2); % 11    loss rate (1/day) at 20 deg C
lake_params{55 -5} = x(3); % 12    specific growth rate (1/day) at 20 deg C
lake_params{58 -5} = x(4); % 15    Half saturation growth P level (mg/m3)
lake_params{61 -5} = x(5); % 18    Settling velocity for Chl2 a (m day-1)
lake_params{62 -5} = x(6);  % 19    Loss rate (1/day) at 20 deg C
lake_params{63 -5} = x(7);  % 20    Specific growth rate (1/day) at 20 deg C
lake_params{64 -5} = x(8);  % 21    Half saturation growth P level (mg/m3)
lake_params{51 -5} = x(9);  % % 8  settling velocity for S (m day-1)
sediment_params{22} = x(10);  % 38 R16 sorption of P on Fe k
sediment_params{34} = x(11);  %    accel

lake_params{21 -5} = x(12); % 16    scaling factor for inflow volume (-)
lake_params{25 -5} = x(13); % 20    scaling factor for inflow concentration of total P (-)
lake_params{26 -5} = x(14); % 21    scaling factor for inflow concentration of diss. organic P (-)
lake_params{39 -5} = x(15); % 34    Scaling factor for inflow concentration of Fe3 (-)



run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID
m_start=[2000, 1, 1]; %
m_stop=[2011, 12, 31]; %

disp(datetime('now'));

[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, is_save_results); % runs the model and outputs obs and sim

zinx=find(MyLake_results.basin1.z<4);
TP_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:)+MyLake_results.basin1.concentrations.PP(zinx,:) + MyLake_results.basin1.concentrations.DOP(zinx,:))', 2);
Chl_mod = mean((MyLake_results.basin1.concentrations.Chl(zinx,:)+MyLake_results.basin1.concentrations.C(zinx,:))', 2);
P_mod = mean((MyLake_results.basin1.concentrations.P(zinx,:))', 2);
PP_mod = mean((MyLake_results.basin1.concentrations.PP(zinx,:))', 2);

load 'obs/store_obs/TOTP.dat' % measured
% load 'obs/store_obs/Cha.dat' % measured
load 'obs/store_obs/Cha_aquaM_march_2017.dat' % measured
load 'obs/store_obs/PO4.dat' % measured
load 'obs/store_obs/Part.dat' % measured


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, TOTP(:,1)));
r_TOTP = RMSE(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha_aquaM_march_2017(:,1)));
r_Chl = RMSE(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
r_PO4 = RMSE(P_mod(loc_sim, 1), PO4(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
r_PP = RMSE(PP_mod(loc_sim, 1), Part(loc_obs, 2));


x'
res = sum([r_TOTP, 2*r_Chl, r_PO4, r_PP])


function r = RMSE(y, yhat)
    r = sqrt(mean((y-yhat).^2));

function [c,ceq] = nonlcon(x)
c = [-x(1); -x(2); -x(3); -x(4); -x(5); -x(6); -x(7); -x(8); -x(9); -x(10); -x(11)];
ceq = [];
