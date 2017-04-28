%% load_params: load parameters for MyLake and Sediment
function [lake_params, sediment_params] = load_params()

lake_params = {
    % PhysPar
    0.5, 'dz',                 % 1
    0.0322, 'Kz_K1',           % 2     open water diffusion parameter (-)
    0.000898, 'Kz_K1_ice',     % 3     under ice diffusion parameter (-)
    2.21E-06, 'Kz_N0',         % 4     min. stability frequency (s-2)
    0.3627, 'C_shelter',       % 5     wind shelter parameter (-)
    59.42, 'lat',              % 6     latitude (decimal degrees)
    10.83, 'lon',              % 7     longitude (decimal degrees)
    0.3, 'alb_melt_ice',       % 8     albedo of melting ice (-)
    0.77, 'alb_melt_snow',     % 9     albedo of melting snow (-)
    2.00E-04, 'PAR_sat',       % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    0.45, 'f_par',             % 11    Fraction of PAR in incoming solar radiation (-)
    0.015, 'beta_chl',         % 12    Optical cross_section of chlorophyll (m2 mg-1)
    5, 'lambda_i',             % 13    PAR light attenuation coefficient for ice (m-1)
    15, 'lambda_s',            % 14    PAR light attenuation coefficient for snow (m-1)
    0.36, 'F_sed_sld',         % 15    volume fraction of solids in sediment (= 1-porosity)
    1, 'I_scV',                % 16    scaling factor for inflow volume (-)
    1, 'I_scT',                % 17    scaling coefficient for inflow temperature (-)
    1, 'I_scC',                % 18    scaling factor for inflow concentration of C (-)
    1, 'I_scS',                % 19    scaling factor for inflow concentration of S (-)
    1, 'I_scTP',             % 20    scaling factor for inflow concentration of total P (-)
    1, 'I_scDOP',              % 21    scaling factor for inflow concentration of diss. organic P (-)
    1, 'I_scChl',              % 22    scaling factor for inflow concentration of Chl a (-)
    1, 'I_scDOC',              % 23    scaling factor for inflow concentration of DOC  (-)
    1, 'I_scPOC',              % 23    scaling factor for inflow concentration of POC  (-)

    % BioPar
    2.5, 'swa_b0',             % 1     non-PAR light attenuation coeff. (m-1)
    1, 'swa_b1',               % 2     PAR light attenuation coeff. (m-1)
    3.30E-07, 'S_res_epi',     % 3     Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
    3.30E-08, 'S_res_hypo',    % 4     Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
    0.03, 'H_sed',             % 5     height of active sediment layer (m, wet mass)
    15, 'Psat_L',              % 6     NOT USED: Half saturation parameter for Langmuir isotherm
    30, 'Fmax_L',              % 7     NOT USED: Scaling parameter for Langmuir isotherm !!!!!!!!!!!!
    0.05, 'w_s',               % 8     settling velocity for S (m day-1)
    0.01, 'w_chl',             % 9     settling velocity for Chl a (m day-1)
    1, 'Y_cp',                 % 10    yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)   1/55*112/1 = 1
    0.2, 'm_twty',             % 11    loss rate (1/day) at 20 deg C
    1.5, 'g_twty',             % 12    specific growth rate (1/day) at 20 deg C
    2.00E-04, 'k_twty',        % 13    NOT USED: specific Chl a to P transformation rate (1/day) at 20 deg C
    0, 'dop_twty',             % 14    NOT USED: specific DOP to P transformation rate (day-1) at 20 deg C
    0.2, 'P_half',             % 15    Half saturation growth P level (mg/m3)
    3.00E-05, 'PAR_sat_2',     % 16    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    0.015, 'beta_chl_2',       % 17    Optical cross_section of chlorophyll (m2 mg-1)
    0.01, 'w_chl_2',           % 18    Settling velocity for Chl a (m day-1)
    0.2, 'm_twty_2',           % 19    Loss rate (1/day) at 20 deg C
    1.5, 'g_twty_2',           % 20    Specific growth rate (1/day) at 20 deg C
    0.2, 'P_half_2',           % 21    Half saturation growth P level (mg/m3)
    0.01, 'oc_DOC',            % 22    Optical cross-section of DOC (m2/mg DOC)
    0.1, 'qy_DOC',             % 23    Quantum yield (mg DOC degraded/mol quanta)
    0.1, 'k_BOD',              % 24    NOT USED: Organic decomposition rate (1/d)
    500, 'k_SOD',              % 25    NOT USED: Sedimentary oxygen demand (mg m-2 d-1)
    1.047, 'theta_bod',        % 26    NOT USED: Temperature adjustment coefficient for BOD, T ? 10 °C
    1.13, 'theta_bod_ice',     % 27    NOT USED: Temperature adjustment coefficient for BOD, T < 10 °C
    1, 'theta_sod',            % 28    NOT USED: Temperature adjustment coefficient for SOD, T ? 10 °C
    1, 'theta_sod_ice',        % 29    NOT USED: Temperature adjustment coefficient for SOD, T < 10 °C
    4, 'BOD_temp_switch',      % 30    NOT USED: Threshold for bod or bod_ice °C
    5.2, 'pH',                 % 31    Lake water pH
    2, 'Mass_Ratio_C_Chl',     % 32    NOT USED: Fixed empirical ratio C:Chl (mass/mass)
    100, 'I_scDIC',            % 33    Scaling factor for inflow concentration of DOC  (-)
    0.25, 'SS_C',              % 34    Carbon fraction in H_netsed_catch
    1.95, 'density_org_H_nc',  % 35    Density of organic fraction in H_netsed_catch [g cm-3]
    2.65, 'density_inorg_H_nc',% 36    Density of inorganic fraction in H_netsed_catch [g cm-3]
    1, 'I_scO',                % 37    Scaling factor for inflow concentration
    2, 'Q10_wc',               % 38    Q10 for reactions of respiration
    1, 'wc_factor',            % 39    Scaling factor for rates in WC
    4.8497, 'T_ref_wc'};       % 40    Reference Temperature for rates

sediment_params = {
    1,     'k_OM1';  % 1
    0.1,  'k_OM2';  % 0.01
    10,   'k_DOM1';  % 0.01
    1,     'k_DOM2';  % 0.01
    0.008, 'Km_O2';     % Canavan, R. W (2006)
    0.01,  'Km_NO3';    % Canavan, R. W (2006)
    0.2,   'Km_Fe(OH)3';  % Canavan, R. W (2006)
    0.2,   'Km_FeOOH';
    0.1,   'Km_SO4';   % Canavan, R. W (2006
    0.001, 'Km_oxao';
    0.1,   'Km_amao';
    0.008, 'Kin_O2';     % the same as Km
    0.01,  'Kin_NO3';    % the same as Km
    0.2,   'Kin_FeOH3';  % the same as Km
    0.2,   'Kin_FeOOH';  % the same as Km
    20,    'k_amox';     % Canavan, R. W (2006)
    5000,  'k_Feox';     % Canavan, R. W (2006)
    0.1,   'k_Sdis';
    2500,  'k_Spre';
    3.17,  'k_FeS2pre';
    0.1,   'k_alum';
    % 1.35,  'k_pdesorb_a';
    300,  'k_pdesorb_a';
    % 1.35,  'k_pdesorb_b';
    300,  'k_pdesorb_b';
    6500,  'k_rhom';
    0.1,   'k_tS_Fe';
    9600,  'Ks_FeS';    % Canavan, R. W (2006)
    0.001, 'k_Fe_dis';  % Canavan, R. W (2006)
    1.5e-3,'k_Fe_pre';  % Canavan, R. W (2006)
    0.37,  'k_apa';
    3e-6,  'kapa';
    0.3134,'k_oms';
    1000,  'k_tsox';   % Canavan, R. W (2006)
    3.3e-3,'k_FeSpre'; % Canavan, R. W (2006)
    30,    'accel';

    1e-6,   'f_pfe';
    1.35,   'k_pdesorb_c';

    % Added porosity modeling parameters:
    0.98,   'fi_in';
    0.85,   'fi_f';
    0.5,    'X_b';
    1,      'tortuosity';

    0.1,    'w';
    256,    'n';
    30,     'depth';
    1,   'F';
    14.4,   'alfa0';

    % OM composition
    106,    'Cx1';
    16,     'Ny1';
    1,      'Pz1';
    200,    'Cx2';
    20,     'Ny2';
    0,      'Pz2';

    0.0001,  'ts';
    };
