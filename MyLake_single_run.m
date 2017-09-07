for i=1:1000
tic
disp('Started at:')
disp(datetime('now'));


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

is_metrics = true; % print metrics in the end

m_start=[2005, 1, 1]; %
m_stop=[2005, 12, 31]; %

save_initial_conditions = false; % save final concentrations as initial for the next run
file_name = 'IO/ecomac_2_res.mat'

[lake_params, sediment_params] = load_params();

% sediment_params{62}  = 14.4/2; % 62 alfa0 bioirrigation
% sediment_params{74}  = 1; % 74 pH algorithm
% lake_params{19}  = 1; % POC scaling
% lake_params{34}  = 50; % Fe3 scaling

% Niva results:
% lake_params{47} = 0.0500000000000000; % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = 0.110668940917121; % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = 1; % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = 0.638922175649520; % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = 0.0593885907392169; % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = 0.114773582911439;  % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = 1.09846108737339;  % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = 1.20000000000000;  % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = 0.0539465787121992;  % % 46  settling velocity for S (m day-1)
% lake_params{23} = 0.521896084255766; % 23    scaling factor for inflow concentration of DOC  (-)
% lake_params{19} = 9.02485516751897; % 19    scaling factor for inflow concentration of POC (-)
% lake_params{34} = 52.0521098578564; % 34    Scaling factor for inflow concentration of Fe3 (-)
% lake_params{10} = 2.45705383730120e-05; % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = 4.73963587003187e-05; % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = 0.0450000000000000; % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = 0.0296430648972882; % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{23} = 62563.6503110866;  % 38 R16 sorption of P on Fe k
% sediment_params{52} = 65.1237288191057;  %    accel
% lake_params{24} = 0.484215194118384;  % 24    scaling factor for inflow concentration of POP (-)
% lake_params{19} = 0;  %    % 19    scaling factor for inflow concentration of POC (-)
% lake_params{34} = 0.0382871701392115;  %    % 34    Scaling factor for inflow concentration of Fe3 (-)


% ecomac-2 results
lake_params{47} = 0.164426980537189; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{49} = 0.254982092942300; % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
lake_params{50} = 1.25000000000000; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{53} = 0.502268862228062; % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
lake_params{56} = 0.217584383024992; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} = 0.100000000000000; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
lake_params{58} = 1.42096457634565; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} = 0.342624466538650; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
lake_params{46} = 0.120588686751271; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{10} = 1.27800593870456e-05; % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} = 1.00000000000000e-05; % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{12} = 0.0450000000000000; % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{55} = 0.0450000000000000; % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
sediment_params{52} = 32.5; % 65.1237e+000   %    accel
lake_params{24} = 1.06941658550207; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)













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
%
