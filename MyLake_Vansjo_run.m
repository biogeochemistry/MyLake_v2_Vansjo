% for i = 1:1000

tic
disp('Started at:')
disp(datetime('now'));


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

is_metrics = true; % print metrics in the end

m_start=[2000, 1, 1]; %
m_stop=[2013, 12, 31]; %
% big_results = cell(1,no_runs);  % collects the results
% big_inputs = cell(1,no_runs);   % collects the inputs
save_initial_conditions = false; % save final concentrations as initial for the next run
file_name = 'IO/test.mat'

[lake_params, sediment_params] = load_params();

% NIVA chl run err=333;
lake_params{47} =  0.0885; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{49} =  0.2553; % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
lake_params{50} =  1.0226; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{53} =  1.6529; % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
lake_params{56} =  0.0795; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} =  0.1205; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
lake_params{58} =  1.1221; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} =  1.0316; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
lake_params{46} =  0.1372; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{10} =  0.0001;  % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} =  0.0001;  % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{12} =  0.0376;  % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{55} =  0.0345;  % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)


% try
run_ID = 0;
clim_ID = 0;
[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim





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
    % load 'obs/store_obs/Cha.dat' % measured
    load 'obs/store_obs/Cha_aquaM_march_2017.dat' % measured
    load 'obs/store_obs/PO4.dat' % measured
    load 'obs/store_obs/Part.dat' % measured


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, TOTP(:,1)));
    rmsd_TOTP = sqrt(mean((TP_mod(loc_sim, 1)-TOTP(loc_obs, 2)).^2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha_aquaM_march_2017(:,1)));
    rmsd_Chl = sqrt(mean((Chl_mod(loc_sim, 1)-Cha_aquaM_march_2017(loc_obs, 2)).^2));


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