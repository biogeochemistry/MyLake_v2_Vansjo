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
m_start=[2002, 1, 1]; %
m_stop=[2011, 12, 31]; %

try
    [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID); % runs the model and outputs obs and sim


    TP_mod = mean((MyLake_results.Pzt(zinx,:)+MyLake_results.PPzt(zinx,:)+MyLake_results.DOPzt(zinx,:)+MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:))', 2);
    Chl_mod = mean((MyLake_results.Chlzt(zinx,:)+MyLake_results.Czt(zinx,:))', 2);
    Pzt_mod = mean((MyLake_results.Pzt(zinx,:))', 2);
    PPzt_mod = mean((MyLake_results.PPzt(zinx,:))', 2);

    load 'obs/store_obs/TOTP.dat' % these are just C&P of vanem ...
    load 'obs/store_obs/Cha.dat' % these are just C&P of vanem ...
    load 'obs/store_obs/PO4.dat' % these are just C&P of vanem ...
    load 'obs/store_obs/Part.dat' % these are just C&P of vanem ...


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, TOTP(:,1)));
    c_TOTP = corrcoef(TP_mod(loc_sim, 1), TOTP(loc_obs, 2), 'alpha',0.05);


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, Cha(:,1)));
    c_Chl = corrcoef(Chl_mod(loc_sim, 1), Cha(loc_obs, 2), 'alpha',0.05);


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, PO4(:,1)));
    c_PO4 = corrcoef(Pzt_mod(loc_sim, 1), PO4(loc_obs, 2), 'alpha',0.05);


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.days, Part(:,1)));
    c_PP = corrcoef(PPzt_mod(loc_sim, 1), Part(loc_obs, 2), 'alpha',0.05);


    res = - mean([c_TOTP(1,2), c_Chl(1,2), c_PO4(1,2), c_PP(1,2)])
catch ME
    res = 1e8
end
