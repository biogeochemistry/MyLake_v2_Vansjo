for i=1:1000
tic
disp('Started at:')
disp(datetime('now'));


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

is_metrics = true; % print metrics in the end

m_start=[2000, 1, 1]; %
m_stop=[2013, 12, 31]; %

save_initial_conditions = true; % save final concentrations as initial for the next run
file_name = 'IO/test_pH_1.mat'

[lake_params, sediment_params] = load_params();

sediment_params{62}  = 14.4/2; % 62 alfa0 bioirrigation
sediment_params{74}  = 1; % 74 pH algorithm
sediment_params{73}  = 25; % (100 is the minimum) number of time steps during 1 day (fixed time step of MyLake) for chemical and sediment module (the modules should be in sync)
lake_params{19}  = 1; % POC scaling

% Niva initial results effective depth 30:
x = [50.0000e-003,   110.6689e-003,     1.0000e+000,   638.9222e-003,   204.8121e-003,   167.6746e-003,     1.0985e+000,     1.5525e+000,    53.9466e-003,   521.8961e-003,    26.4624e+000,    51.0521e+000,    24.5705e-006,    75.5867e-006,    45.0000e-003,    29.6431e-003,     4.9321e+003,    65.1237e+000,   390.1162e-003,     0.0000e+000,    84.3021e-003]

lake_params{47} = x(1); % 47     settling velocity for Chl1 a (m day-1)
lake_params{49} = x(2); % 49    loss rate (1/day) at 20 deg C
lake_params{50} = x(3); % 50    specific growth rate (1/day) at 20 deg C
lake_params{53} = x(4); % 53    Half saturation growth P level (mg/m3)
lake_params{56} = x(5); % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} = x(6);  % 57    Loss rate (1/day) at 20 deg C
lake_params{58} = x(7);  % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} = x(8);  % 59    Half saturation growth P level (mg/m3)
lake_params{46} = x(9);  % % 46  settling velocity for S (m day-1)
lake_params{23} = x(10); % 23    scaling factor for inflow concentration of DOC  (-)
lake_params{19} = x(11); % 19    scaling factor for inflow concentration of POC (-)
lake_params{34} = x(12); % 34    Scaling factor for inflow concentration of Fe3 (-)
lake_params{10} = x(13); % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} = x(14); % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{12} = x(15); % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{55} = x(16); % 17    Optical cross_section of chlorophyll (m2 mg-1)
sediment_params{23} = x(17);  % 38 R16 sorption of P on Fe k
sediment_params{52} = x(18);  %    accel
lake_params{24} = x(19);  % 24    scaling factor for inflow concentration of POP (-)
lake_params{19} = x(20);  %    % 19    scaling factor for inflow concentration of POC (-)
lake_params{34} = x(21);  %    % 34    Scaling factor for inflow concentration of Fe3 (-)



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
end
