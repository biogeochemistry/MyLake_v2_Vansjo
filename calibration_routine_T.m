function x = calibration_routine()
tic
format shortEng
format compact
% parpool
% gaoptions = optimoptions('ga','UseParallel',true);

% x = [0.032; 0.3627];
% lb = x0*0.1;
% ub = x0*10;

x(1) = 44.2401e-003;; %, 'Kz_K1',           % 2     open water diffusion parameter (-)
x(2) = 0.000898; %, 'Kz_K1_ice',     % 3     under ice diffusion parameter (-)
x(3) = 7E-05; %, 'Kz_N0',         % 4     min. stability frequency (s-2)
x(4) = 0.5; %, 'C_shelter',       % 5     wind shelter parameter (-)
x(5) = 0.3; %, 'alb_melt_ice',       % 8     albedo of melting ice (-)
x(6) = 0.77; %, 'alb_melt_snow',     % 9     albedo of melting snow (-)

lb = zeros(size(x));
ub = ones(size(x));

% lb = [0.01; 0.1];
% ub = [0.1; 0.5];


fcns = {@gaplotscorediversity, @gaplotstopping, @gaplotgenealogy, @gaplotscores, @gaplotdistance, @gaplotselection, @gaplotmaxconstr, @gaplotbestf, @gaplotbestindiv, @gaplotexpectation, @gaplotrange, @gaplotpareto, @gaplotparetodistance, @gaplotrankhist, @gaplotspread};

population_size = 72;  % Populations size for each generation of the genetic algorithm
max_generations = 7;  % How many generations to run the genetic algorithm for
parallelize     = true; % 15 generation takes 12 hours on 24 cores

% options = gaoptimset('Display','iter','UseParallel', true, 'TolFun', 1e-2, 'PlotFcns', fcns);
options = optimoptions('ga', 'MaxGenerations', max_generations, 'PopulationSize', population_size, 'UseParallel', parallelize);

x = ga(@opt_fun, length(x) ,[],[],[],[],lb,ub, @nonlcon, options)

%% opt_fun: function which we are going to minimize
function [res] = opt_fun(x)

[lake_params, sediment_params] = load_params();


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files
is_save_results = false; % Do not save profiles and initial concentrations



lake_params{2} = x(1); % 2     open water diffusion parameter (-)
lake_params{3} = x(2); % 3     under ice diffusion parameter (-)
lake_params{4} = x(3); % 4     min. stability frequency (s-2)
lake_params{5} = x(4); % 5     wind shelter parameter (-)
lake_params{8} = x(5); % 8     albedo of melting ice (-)
lake_params{9} = x(6); % 9     albedo of melting snow (-)




run_ID = 'Vansjo_Hist_M0' ; %  CALIBRATION RUN
clim_ID = run_ID
m_start=[2004, 1, 1]; %
m_stop=[2013, 12, 31]; %

disp(datetime('now'));
[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, is_save_results); % runs the model and outputs obs and sim


load('Postproc_code/Vansjo/VAN1_data_2017_02_28_10_55.mat')

depths = [5;10;15;20;25;30;35;40];
r_Temp = 0;


for i=1:size(depths,1)
    d = depths(i);
    zinx=find(MyLake_results.basin1.z == d);
    T_measured = res.T(res.depth1 == d);
    day_measured = res.date(res.depth1 == d);
    day_measured = day_measured(~isnan(T_measured));
    T_measured = T_measured(~isnan(T_measured));

    Temp_mod = MyLake_results.basin1.T(zinx,:)';
    [T_date,loc_sim, loc_obs] = intersect(MyLake_results.basin1.days, day_measured);

    r_Temp = r_Temp + RMSE(Temp_mod(loc_sim, 1), T_measured(loc_obs, 1));
end

x'
res = r_Temp



function r = RMSE(y, yhat)
    r = sqrt(mean((y-yhat).^2));

function [c,ceq] = nonlcon(x)
c = [-x(1); -x(2)]; % -x(3); -x(4); -x(5); -x(6); -x(7); -x(8); -x(9); -x(10); -x(11)];
ceq = [];
