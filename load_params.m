%% load_params: load parameters for MyLake and Sediment
function [lake_params, sediment_params] = load_params()

lake_params = {
    % PhysPar
    0.5, 'dz',                 % 1
    0.0164, 'Kz_K1',           % 2     open water diffusion parameter (-)
    0.000898, 'Kz_K1_ice',     % 3     under ice diffusion parameter (-)
    7E-05, 'Kz_N0',            % 4     min. stability frequency (s-2)
    0.74, 'C_shelter',         % 5     wind shelter parameter (-)
    59.40, 'lat',              % 6     latitude (decimal degrees)
    10.80, 'lon',              % 7     longitude (decimal degrees)
    0.3, 'alb_melt_ice',       % 8     albedo of melting ice (-)
    0.77, 'alb_melt_snow',     % 9     albedo of melting snow (-)
    2.00E-04, 'PAR_sat',       % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    0.45, 'f_par',             % 11    Fraction of PAR in incoming solar radiation (-)
    0.015, 'beta_chl',         % 12    Optical cross_section of chlorophyll (m2 mg-1)
    5, 'lambda_i',             % 13    PAR light attenuation coefficient for ice (m-1)
    15, 'lambda_s',            % 14    PAR light attenuation coefficient for snow (m-1)
    0.36, 'F_sed_sld',         % 15    volume fraction of solids in sediment (= 1-porosity)
    1, 'I_scV',                % 16    scaling factor for inflow volume (-)
    0, 'I_scT',                % 17    adjusting delta for inflow temperature (-)
    1, 'I_scC',                % 18    scaling factor for inflow concentration of C (-)
    1, 'I_scPOC',              % 19    scaling factor for inflow concentration of POC (-)
    1, 'I_scTP',               % 20    scaling factor for inflow concentration of total P (-)
    1, 'I_scDOP',              % 21    scaling factor for inflow concentration of diss. organic P (-)
    1, 'I_scChl',              % 22    scaling factor for inflow concentration of Chl a (-)
    1, 'I_scDOC',              % 23    scaling factor for inflow concentration of DOC  (-)
    1, 'I_scPOP',              % 24    scaling factor for inflow concentration of POP  (-)
    1, 'I_scO',                % 25    Scaling factor for inflow concentration of O2 (-)
    1 , 'I_scDIC',             % 26    Scaling factor for inflow concentration of DOC  (-)
    1,  'I_scNO3',             % 27    Scaling factor for inflow concentration of NO3 (-)
    1,  'I_scNH4',             % 28    Scaling factor for inflow concentration of NH4 (-)
    1,  'I_scSO4',             % 29    Scaling factor for inflow concentration of SO4 (-)
    1,  'I_scFe2',             % 30    Scaling factor for inflow concentration of Fe2 (-)
    1,  'I_scCa2',             % 31    Scaling factor for inflow concentration of Ca2 (-)
    1,  'I_scpH',              % 32    Scaling factor for inflow concentration of pH (-)
    1,  'I_scCH4',             % 33    Scaling factor for inflow concentration of CH4 (-)
    1,  'I_scFe3',             % 34    Scaling factor for inflow concentration of Fe3 (-)
    1,  'I_scAl3',             % 35    Scaling factor for inflow concentration of Al3 (-)
    1,  'I_scSiO4',            % 36    Scaling factor for inflow concentration of SiO4 (-)
    1,  'I_scSiO2',            % 37    Scaling factor for inflow concentration of SiO2 (-)
    1,  'I_scCH4g',            % 38    Scaling factor for inflow concentration of CH4g (-)
    2.5, 'swa_b0',             % 39     non-PAR light attenuation coeff. (m-1)
    1.05, 'swa_b1',            % 40     PAR light attenuation coeff. (m-1)
    3.30E-07, 'S_res_epi',     % 41     Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
    3.30E-08, 'S_res_hypo',    % 42     Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
    0.03, 'H_sed',             % 43     height of active sediment layer (m, wet mass)
    15, 'Psat_L',              % 44     NOTE: NOT USED: Half saturation parameter for Langmuir isotherm
    30, 'Fmax_L',              % 45     NOTE: NOT USED: Scaling parameter for Langmuir isotherm !!!!!!!!!!!!
    0.05, 'w_s',               % 46     settling velocity for S (m day-1)
    0.01, 'w_chl',             % 47     settling velocity for Chl a (m day-1)
    1, 'Y_cp',                 % 48     NOTE: NOT USED:  yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)   1/55*112/1 = 1
    0.2, 'm_twty',             % 49    loss rate (1/day) at 20 deg C
    1.5, 'g_twty',             % 50    specific growth rate (1/day) at 20 deg C
    2.00E-04, 'k_twty',        % 51    NOTE: NOT USED: specific Chl a to P transformation rate (1/day) at 20 deg C
    0, 'dop_twty',             % 52    NOTE: NOT USED: specific DOP to P transformation rate (day-1) at 20 deg C
    0.2, 'P_half',             % 53    Half saturation growth P level (mg/m3)
    3.00E-05, 'PAR_sat_2',     % 54    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    0.015, 'beta_chl_2',       % 55    Optical cross_section of chlorophyll (m2 mg-1)
    0.01, 'w_chl_2',           % 56    Settling velocity for Chl a (m day-1)
    0.2, 'm_twty_2',           % 57    Loss rate (1/day) at 20 deg C
    1.5, 'g_twty_2',           % 58    Specific growth rate (1/day) at 20 deg C
    0.2, 'P_half_2',           % 59    Half saturation growth P level (mg/m3)
    0.01, 'oc_DOC',            % 60    Optical cross-section of DOC (m2/mg DOC)
    0.1, 'qy_DOC',             % 61    Quantum yield (mg DOC degraded/mol quanta)
    0.1, 'k_BOD',              % 62    NOTE: NOT USED: Organic decomposition rate (1/d)
    5, 'w_CH4',                % 63    Methane gas rising velocity (m/d)
    1.047, 'theta_bod',        % 64    NOTE: NOT USED: Temperature adjustment coefficient for BOD, T ? 10 °C
    1.13, 'theta_bod_ice',     % 65    NOTE: NOT USED: Temperature adjustment coefficient for BOD, T < 10 °C
    1, 'theta_sod',            % 66    NOTE: NOT USED: Temperature adjustment coefficient for SOD, T ? 10 °C
    1, 'theta_sod_ice',        % 67    NOTE: NOT USED: Temperature adjustment coefficient for SOD, T < 10 °C
    4, 'BOD_temp_switch',      % 68    NOTE: NOT USED: Threshold for bod or bod_ice °C
    7.5, 'pH',                 % 69    Lake water pH
    2, 'Q10_wc',               % 70    Q10 for reactions of respiration
    1, 'wc_factor',            % 71    Scaling factor for rates in WC
    4.8497, 'T_ref_wc'};       % 72    Reference Temperature for rates

sediment_params = {
    10,    'k_Chl';        % 1       % 1
    1,     'k_POP';        % 2       % 1
    1      'k_POC';        % 3       % 0.01
    1,     'k_DOP';        % 4       % 1
    1,     'k_DOC';        % 5       % 1
    0.008, 'Km_O2';        % 6       % Canavan, R. W (2006) rho=2
    0.01,  'Km_NO3';       % 7       % Canavan, R. W (2006) rho=2
    0.2*2,  'Km_Fe(OH)3';   % 8       % Canavan, R. W (2006) rho=2
    0.2*2,  'Km_FeOOH';     % 9       % Canavan, R. W (2006) rho=2
    0.1,  'Km_SO4';       % 10       % Canavan, R. W (2006 rho=2
    0.001,'Km_oxao';      % 11       % the same as Km rho=2
    0.1,  'Km_amao';      % 12       % the same as Km rho=2
    0.008, 'Kin_O2';       % 13       % the same as Km rho=2
    0.01,  'Kin_NO3';      % 14       % the same as Km rho=2
    0.2*2,   'Kin_FeOH3';    % 15       % the same as Km rho=2
    0.2*2,   'Kin_FeOOH';    % 16       % the same as Km rho=2
    20,    'k_amox';       % 17       % Canavan, R. W (2006)
    50e3,  'k_Feox';       % 18       % Canavan, R. W (2006)
    0.1,   'k_Sdis';       % 19       %
    2500,  'k_Spre';       % 20       %
    3.17,  'k_FeS2pre';    % 21
    0.1,   'k_alum';       % 22
    1.35,  'k_pdesorb_a';; % 23
    1.35,  'k_pdesorb_b';  % 24
    6500,  'k_rhom';       % 25
    8,   'k_tS_Fe';        % 26      % Cappellen (1996) in Canavan, R. W (2006) the reaction is different
    9600,  'Ks_FeS';       % 27      % Canavan, R. W (2006)
    0.001, 'k_Fe_dis';     % 28      % Canavan, R. W (2006)
    1500*2,'k_Fe_pre';     % 29      % Canavan, R. W (2006)
    0.37,  'k_apa';        % 30
    3e-6,  'kapa';         % 31
    0.3134,'k_oms';        % 32
    1e6,   'k_tsox';       % 33     % Canavan, R. W (2006)
    0.3/2, 'k_FeSpre';     % 34     % from "Non-steady state diagenesis of organic and inorganic sulfur in lake sediments Raoul-Marie Couture, Rachele Fischer b, Philippe Van Cappellen b, Charles Gobeil c
    1e7,   'k_ch4_o2'      % 35     % Canavan, R. W (2006)
    1e-1,  'k_ch4_so4'     % 36     % Canavan, R. W (2006)
    2.66,  'CH4_solubility' % 37     % Canavan, R. W (2006)
    1e3,   'k_ch4_dis'             % 38
    5,     'CH4_rising_vel' % 39     % Canavan, R. W (2006)
    30,    'accel';        % 40
    1e-6,   'f_pfe';       % 41
    1.35,   'k_pdesorb_c'; % 42
    0.98,   'fi_in';       % 43
    0.85,   'fi_f';        % 44
    0.5,    'X_b';         % 45
    1,      'tortuosity';  % 46
    0.1,    'w';           % 47
    256,    'n';           % 48
    30,     'depth';       % 49
    5,      'w_CH4';       % 50     rising velocity of methane gas in the sediment
    14.4,   'alfa0';       % 51
    106,    'Cx1';         % 52           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    16,     'Ny1';         % 53           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Pz1';         % 54           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    200,    'Cx2';         % 55           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    20,     'Ny2';         % 56           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Pz2';         % 57           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Cx3';         % 58           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    0.1,    'Ny3';         % 59           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    0,      'Pz3';         % 60           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    10,      'n_ts';       % 61           % number of time steps during 1 day (fixed time step of MyLake) for chemical and sediment module (the modules should be in sync)
    };
