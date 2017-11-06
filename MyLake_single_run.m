for i=1:1000
tic
disp('Started at:')
disp(datetime('now'));


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

is_metrics = true; % print metrics in the end

m_start=[1995, 1, 1]; %
m_stop=[2013, 12, 31]; %

save_initial_conditions = false; % save final concentrations as initial for the next run


[lake_params, sediment_params] = load_params();

% % Niva results RMSD = 130 =======================================
% file_name = 'IO/calibrated.mat'
% lake_params{47} = 58.3842e-003; % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
% lake_params{49} = 128.2949e-003; % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
% lake_params{50} = 1.4988e+000; % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
% lake_params{53} = 1.6945e+000; % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
% lake_params{56} = 208.3324e-003; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
% lake_params{57} = 201.6135e-003; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
% lake_params{58} = 1.2687e+000; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
% lake_params{59} = 1.6142e+000; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
% lake_params{46} = 31.3665e-003; % 53.9466e-003   % % 46  settling velocity for S (m day-1)
% lake_params{10} = 14.4699e-006; % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{54} = 30.5827e-006; % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{12} = 37.9560e-003; % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
% lake_params{55} = 34.7141e-003; % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
% sediment_params{52} = 21.5114e+000; % 65.1237e+000   %    accel
% lake_params{24} = 373.1228e-003; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)

% % Trials:
% lake_params{24} = 1; % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
% =====================================================================================================================


% Latest Niva calibration with RMSD 133: % ====================================================
file_name = 'IO/calib_17_11_01_Niva_calib_rmsd_3x_1x_133.mat'
x = [0.0691969351707719; 0.106384525900796; 1.33423213718069; 1.81948427216433; 0.0992004077859556; 0.211874114480601; 1.08692651868250; 0.662212213886266; 0.0426104938223779; 2.28258392431412e-05; 7.44071863750726e-05; 0.0432478460870351; 0.0232570145845245; 22.0651718079512; 0.735236639414830];
lake_params{47} = x(1); % 50.0000e-003  % 47     settling velocity for Chl1 a (m day-1)
lake_params{49} = x(2); % 110.6689e-003  % 49    loss rate (1/day) at 20 deg C
lake_params{50} = x(3); % 1.0000e+000  % 50    specific growth rate (1/day) at 20 deg C
lake_params{53} = x(4); % 638.9222e-003  % 53    Half saturation growth P level (mg/m3)
lake_params{56} = x(5); % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} = x(6); % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
lake_params{58} = x(7); % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} = x(8); % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
lake_params{46} = x(9); % 53.9466e-003   % % 46  settling velocity for S (m day-1)
lake_params{10} = x(10); % 24.5705e-006  % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{54} = x(11); % 75.5867e-006  % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
lake_params{12} = x(12); % 45.0000e-003  % 12    Optical cross_section of chlorophyll (m2 mg-1)
lake_params{55} = x(13); % 29.6431e-003  % 17    Optical cross_section of chlorophyll (m2 mg-1)
sediment_params{52} = x(14); % 65.1237e+000   %    accel
lake_params{24} = x(15); % 390.1162e-003   % 24    scaling factor for inflow concentration of POP (-)
% % =====================================================================================================================



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
    rmsd_TOTP = rmsd(TP_mod(loc_sim, 1), TOTP(loc_obs, 2))


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Cha_aquaM_march_2017(:,1)));
    rmsd_Chl = rmsd(Chl_mod(loc_sim, 1), Cha_aquaM_march_2017(loc_obs, 2));
    

    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, PO4(:,1)));
    rmsd_PO4 = rmsd(P_mod(loc_sim, 1), PO4(loc_obs, 2));


    [TP_date,loc_sim, loc_obs] = (intersect(MyLake_results.basin1.days, Part(:,1)));
    rmsd_PP = rmsd(POP_mod(loc_sim, 1), Part(loc_obs, 2));

    disp('RMSD 3xRMSE(P)+RMSE(O2):')
    disp(sum([3*rmsd_TOTP, 3*rmsd_Chl, 3*rmsd_PO4, 3*rmsd_PP, rmsd_O2]))
    disp('RMSD = RMSE(P)+RMSE(O2):')
    disp(sum([rmsd_TOTP, rmsd_Chl, rmsd_PO4, rmsd_PP, rmsd_O2]))
end


toc
end
%
