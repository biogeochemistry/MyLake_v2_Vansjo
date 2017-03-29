function x = calibration_routine()
tic
format shortEng
format compact
% parpool
% gaoptions = optimoptions('ga','UseParallel',true);
x0 = [2500; 8000; 0.02; 0.2; 1.5; 0.2; 0.02;  0.2; 1.5; 0.2; 2];
lb = x0*0.5;
ub = x0*1.5;

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
use_INCA = 1; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files


lake_params{36-7} = x(1); % 6     Half saturation parameter for Langmuir isotherm
lake_params{37-7} = x(2); % 7     Scaling parameter for Langmuir isotherm !!!!!!!!!!!!
lake_params{39-7} = x(3); % 9     settling velocity for Chl1 a (m day-1)
lake_params{41-7} = x(4); % 11    loss rate (1/day) at 20 deg C
lake_params{42-7} = x(5); % 12    specific growth rate (1/day) at 20 deg C
lake_params{45-7} = x(6); % 15    Half saturation growth P level (mg/m3)
lake_params{48-7} = x(7); % 18    Settling velocity for Chl2 a (m day-1)
lake_params{49-7} = x(8);  % 19    Loss rate (1/day) at 20 deg C
lake_params{50-7} = x(9);  % 20    Specific growth rate (1/day) at 20 deg C
lake_params{51-7} = x(10);  % 21    Half saturation growth P level (mg/m3)
lake_params{68-7} = x(11);  % 38    Q10 for reactions of respiration


run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID
m_start=[2000, 1, 1]; %
m_stop=[2011, 12, 31]; %

[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID); % runs the model and outputs obs and sim

zinx=find(MyLake_results.zz<4);
TP_mod = mean((MyLake_results.Pzt(zinx,:)+MyLake_results.PPzt(zinx,:)+MyLake_results.DOPzt(zinx,:)+MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:))', 2);
Chl_mod = mean((MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:))', 2);
Pzt_mod = mean((MyLake_results.Pzt(zinx,:))', 2);
PPzt_mod = mean((MyLake_results.PPzt(zinx,:))', 2);

load 'obs/store_obs/TOTP.dat' % measured
% load 'obs/store_obs/Cha.dat' % measured
load 'obs/store_obs/Cha_aquaM_march_2017.dat' % measured
load 'obs/store_obs/PO4.dat' % measured
load 'obs/store_obs/Part.dat' % measured


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, TOTP(:,1)));
c_TOTP = RMSE(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, Cha_aquaM_march_2017(:,1)));
c_Chl = RMSE(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, PO4(:,1)));
c_PO4 = RMSE(Pzt_mod(loc_sim, 1), PO4(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, Part(:,1)));
c_PP = RMSE(PPzt_mod(loc_sim, 1), Part(loc_obs, 2));


x
res = sum([c_TOTP, c_Chl, c_PO4, c_PP])


function r = RMSE(y, yhat)
    r = sqrt(mean((y-yhat).^2));

function [c,ceq] = nonlcon(x)
c = [-x(1); -x(2); -x(3); -x(4); -x(5); -x(6); -x(7); -x(8); -x(9); -x(10); -x(11)];
ceq = [];
