function x = calibration_routine()
tic
format shortEng
format compact
% parpool
% gaoptions = optimoptions('ga','UseParallel',true);
x = [245.0402e-003, 193.8201e-003, 1.3948e+000, 700.8806e-003, 76.1493e-003, 121.8722e-003, 1.3160e+000, 363.7310e-003, 77.1076e-003, 443.2972e+000, 34.9491e+000];
% lb = x0*0.1;
% ub = x0*10;

lb = [0.05; 0.1; 1; 0.2; 0.05; 0.1; 1; 0.2; 0.01; 1; 1]
ub = [0.5; 0.3; 1.5; 2; 0.5; 0.3; 1.5; 2; 1; 1000; 100]


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


lake_params{40 -7} = x(1); % 9     settling velocity for Chl1 a (m day-1)
lake_params{42 -7} = x(2); % 11    loss rate (1/day) at 20 deg C
lake_params{43 -7} = x(3); % 12    specific growth rate (1/day) at 20 deg C
lake_params{46 -7} = x(4); % 15    Half saturation growth P level (mg/m3)
lake_params{49 -7} = x(5); % 18    Settling velocity for Chl2 a (m day-1)
lake_params{50 -7} = x(6);  % 19    Loss rate (1/day) at 20 deg C
lake_params{51 -7} = x(7);  % 20    Specific growth rate (1/day) at 20 deg C
lake_params{52 -7} = x(8);  % 21    Half saturation growth P level (mg/m3)
lake_params{39 -7} = x(9);  %   settling velocity for S (m day-1)
sediment_params{37-7} = x(10);  %    R16 sorption of P on Fe k
sediment_params{34} = x(11);  %    accel



run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID
m_start=[2000, 1, 1]; %
m_stop=[2011, 12, 31]; %


disp(datetime('now'));
[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID); % runs the model and outputs obs and sim

zinx=find(MyLake_results.z<4);
TP_mod = mean((MyLake_results.Pzt(zinx,:)+MyLake_results.PPzt(zinx,:) + MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:)+MyLake_results.DOPzt(zinx,:))', 2);
Chl_mod = mean((MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:))', 2);
Pzt_mod = mean((MyLake_results.Pzt(zinx,:))', 2);
PPzt_mod = mean((MyLake_results.PPzt(zinx,:))', 2);

load 'obs/store_obs/TOTP.dat' % measured
% load 'obs/store_obs/Cha.dat' % measured
load 'obs/store_obs/Cha_aquaM_march_2017.dat' % measured
load 'obs/store_obs/PO4.dat' % measured
load 'obs/store_obs/Part.dat' % measured


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, TOTP(:,1)));
r_TOTP = RMSE(TP_mod(loc_sim, 1), TOTP(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, Cha_aquaM_march_2017(:,1)));
r_Chl = RMSE(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, PO4(:,1)));
r_PO4 = RMSE(Pzt_mod(loc_sim, 1), PO4(loc_obs, 2));


[TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, Part(:,1)));
r_PP = RMSE(PPzt_mod(loc_sim, 1), Part(loc_obs, 2));


x'
res = sum([r_TOTP, 2*r_Chl, r_PO4, r_PP])


function r = RMSE(y, yhat)
    r = sqrt(mean((y-yhat).^2));

function [c,ceq] = nonlcon(x)
c = [-x(1); -x(2); -x(3); -x(4); -x(5); -x(6); -x(7); -x(8); -x(9); -x(10); -x(11)];
ceq = [];
