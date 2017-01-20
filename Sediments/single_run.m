function [ sediment_results ] = single_run( )
% clear all; close all; clc;
%RUN Summary of this function goes here
%   Detailed explanation goes here

global sed_par_file


% Change here:
% ==============================================
T = 365; % how many days to run
% WC params:
Temperature = 8; % Temperature at SWI
z_max   = 12; % depth at SWI
pH_SWI = 5.2; % pH

params_sediment ...
    = { 0,  'k_OM1';
    0,  'k_OM2';
    0.0123,'Km_O2';
    0.01,  'Km_NO3';
    3.92,  'Km_Fe(OH)3';
    2415,  'Km_FeOOH';
    0.0293,'Km_SO4';
    0.001, 'Km_oxao';
    0.1,   'Km_amao';
    0.3292,'Kin_O2';
    0.1,   'Kin_NO3';
    0.1,   'Kin_FeOH3';
    0.1,   'Kin_FeOOH';
    0  'k_NH4ox';
    0, 'k_Feox';
    0,   'k_Sdis';
    0  'k_Spre';
    0,  'k_FeS2pre';
    0,   'k_alum';
    0,  'k_pdesorb_a';
    0,  'k_pdesorb_b';
    0  'k_rhom';
    0,   'k_tS_Fe';
    2510,  'Ks_FeS';
    0, 'k_Fe_dis';
    0,  'k_Fe_pre';
    0,  'k_apa';
    3e-6,  'kapa';
    0,'k_oms';
    0  'k_tsox';
    0, 'k_FeSpre';
    30,    'accel';
    1e-6,   'f_pfe';
    0,   'k_pdesorb_c';

    % Added porosity modeling parameters:
    0.99999,   'fi_in';
    0.99999,   'fi_f';
    0.5,    'X_b';
    1,      'tortuosity';

    1,      'w';
    128,     'n';
    100,     'depth';
    0.26,   'F';
    0,   'alfa0';

    % OM composition
    112,    'Cx1';
    10,     'Ny1';
    1,      'Pz1';
    200,    'Cx2';
    20,     'Ny2';
    1,      'Pz2';

    0.0001,  'ts'};
% ===================================================


calibration_k_values = [(1:50)',cell2mat(params_sediment(:,1)) ]; % writing sediment

sed_par_file = tempname;
dlmwrite(sed_par_file, calibration_k_values,'delimiter','\t');



[sediment_concentrations, sediment_params, sediment_matrix_templates, species_sediment]  = sediment_init( pH_SWI, z_max, Temperature );

for i=1:T

    % Change here:
    % ==========================================================
    % You may specify the top Boundary conditions here. they can be varied as f(T):
    sediment_bc = {...
        1,      'Ox_c';
        60,        'OM1_fx';
        30,        'OM2_fx';
        0.5,       'PO4_c';
        10,        'NO3_c';
        14.7,      'FeOH3_fx';
        10,        'SO4_c';
        10,        'Fe2_c';
        0,         'FeOOH_fx';
        0,         'FeS_fx';
        0,         'S0_c';
        0,         'S8_fx';
        0,         'FeS2_fx';
        0,         'AlOH3_fx';
        0,         'PO4adsa_fx';
        0,         'PO4adsb_fx';
        10,        'Ca2_c';
        0,         'Ca3PO42_fx';
        0,         'OMS_fx';

        % These values cannot be 0:
        0.1,      'H_c';
        0.1,      'OH_c';
        100,      'CO2_c';
        100,      'CO3_c';
        100,      'HCO3_c';
        0.1,      'NH3_c';
        0.1,      'NH4_c';
        0.1,      'HS_c';
        1E-10,    'H2S_c';
        100,      'H2CO3_c';
    };
    % ========================================================



    sediment_bc = containers.Map({sediment_bc{:,2}},{sediment_bc{:,1}});


    % Running sediment module
    [sediment_bioirrigation_fluxes, sediment_SWI_fluxes, sediment_integrated_over_depth_fluxes, sediment_concentrations, z_sediment, R_values_sedimentz] = sediment_v2(...
    sediment_concentrations, sediment_params, sediment_matrix_templates, species_sediment, sediment_bc);


    % Output:
    O2_sediment_zt(:,i) = sediment_concentrations('Oxygen');
    OM_sediment_zt(:,i) = sediment_concentrations('OM1');
    OMb_sediment_zt(:,i) = sediment_concentrations('OM2');
    NO3_sediment_zt(:,i) = sediment_concentrations('NO3');
    FeOH3_sediment_zt(:,i) = sediment_concentrations('FeOH3');
    SO4_sediment_zt(:,i) = sediment_concentrations('SO4');
    NH4_sediment_zt(:,i) = sediment_concentrations('NH4');
    Fe2_sediment_zt(:,i) = sediment_concentrations('Fe2');
    FeOOH_sediment_zt(:,i) = sediment_concentrations('FeOOH');
    H2S_sediment_zt(:,i) = sediment_concentrations('H2S');
    HS_sediment_zt(:,i)  = sediment_concentrations('HS');
    FeS_sediment_zt(:,i) = sediment_concentrations('FeS');
    S0_sediment_zt(:,i)  = sediment_concentrations('S0');
    PO4_sediment_zt(:,i) = sediment_concentrations('PO4');
    S8_sediment_zt(:,i) = sediment_concentrations('S8');
    FeS2_sediment_zt(:,i) = sediment_concentrations('FeS2');
    AlOH3_sediment_zt(:,i) = sediment_concentrations('AlOH3');
    PO4adsa_sediment_zt(:,i) = sediment_concentrations('PO4adsa');
    PO4adsb_sediment_zt(:,i) = sediment_concentrations('PO4adsb');
    Ca2_sediment_zt(:,i) = sediment_concentrations('Ca2');
    Ca3PO42_sediment_zt(:,i) = sediment_concentrations('Ca3PO42');
    OMS_sediment_zt(:,i) = sediment_concentrations('OMS');
    H_sediment_zt(:,i) = sediment_concentrations('H');
    OH_sediment_zt(:,i) = sediment_concentrations('OH');
    CO2_sediment_zt(:,i) = sediment_concentrations('CO2');
    CO3_sediment_zt(:,i) = sediment_concentrations('CO3');
    HCO3_sediment_zt(:,i) = sediment_concentrations('HCO3');
    NH3_sediment_zt(:,i) = sediment_concentrations('NH3');
    H2CO3_sediment_zt(:,i) = sediment_concentrations('H2CO3');
    pH_sediment_zt(:,i) = -log10(H_sediment_zt(:,i)*10^-6);
    O2_flux_sediment_zt(i) = sediment_SWI_fluxes{1};
    OM_flux_sediment_zt(i) = sediment_SWI_fluxes{2};
    OM2_flux_sediment_zt(i) = sediment_SWI_fluxes{3};
    PO4_flux_sediment_zt(i) = sediment_SWI_fluxes{4};
    % R1_sediment_zt(:,i) = R_values_sedimentz{1};
    % R1_int_sediment_zt(:,i) = R_values_sedimentz{2};
    % R2_sediment_zt(:,i) = R_values_sedimentz{3};
    % R2_int_sediment_zt(:,i) = R_values_sedimentz{4};
    % R3_sediment_zt(:,i) = R_values_sedimentz{5};
    % R3_int_sediment_zt(:,i) = R_values_sedimentz{6};
    % R4_sediment_zt(:,i) = R_values_sedimentz{7};
    % R4_int_sediment_zt(:,i) = R_values_sedimentz{8};
    % R5_sediment_zt(:,i) = R_values_sedimentz{9};
    % R5_int_sediment_zt(:,i) = R_values_sedimentz{10};
    O2_Bioirrigation_sedimentz(:,i) = sediment_bioirrigation_fluxes{1};
    PO4_Bioirrigation_sedimentz(:,i) = sediment_bioirrigation_fluxes{2};
    NO3_flux_sediment_zt(i) = sediment_SWI_fluxes{5};
    FeOH3_flux_sediment_zt(i) = sediment_SWI_fluxes{6};
    % R6_sediment_zt(:,i) = R_values_sedimentz{11};
    % R6_int_sediment_zt(:,i) = R_values_sedimentz{12};

end

% R_values_sediment_zt = {

%       R1_sediment_zt,         'R1';
%       R1_int_sediment_zt,     'R1 integrated';
%       R2_sediment_zt,         'R2';
%       R2_int_sediment_zt,     'R2 integrated';
%       R3_sediment_zt,         'R3';
%       R3_int_sediment_zt,     'R3 integrated';
%       R4_sediment_zt,         'R4';
%       R4_int_sediment_zt,     'R4 integrated';
%       R5_sediment_zt,         'R5';
%       R5_int_sediment_zt,     'R5 integrated';
%       R6_sediment_zt,         'R6';
%       R6_int_sediment_zt,     'R6 integrated';
% };

Bioirrigation_sediment_zt = { O2_Bioirrigation_sedimentz,  'Oxygen';
                              PO4_Bioirrigation_sedimentz, 'PO4'};


sediment_results = {O2_sediment_zt,     'Oxygen (aq)';
                   FeOH3_sediment_zt,   'Iron hydroxide pool 1 Fe(OH)3 (s)';
                   FeOOH_sediment_zt,   'Iron Hydroxide pool 2 FeOOH (s)';
                   SO4_sediment_zt,     'Sulfate SO4(2-) (aq)';
                   Fe2_sediment_zt,     'Iron Fe(2+) (aq)';
                   H2S_sediment_zt,     'Sulfide H2S (aq)';
                   HS_sediment_zt,      'Sulfide HS(-) (aq)';
                   FeS_sediment_zt,     'Iron Sulfide FeS (s)';
                   OM_sediment_zt,      'Organic Matter pool 1 OMa (s)';
                   OMb_sediment_zt,     'Organic Matter pool 2 OMb (s)';
                   OMS_sediment_zt,     'Sulfured Organic Matter (s)';
                   AlOH3_sediment_zt,   'Aluminum oxide Al(OH)3 (s)';
                   S0_sediment_zt,      'Elemental sulfur S(0) (aq)';
                   S8_sediment_zt,      'Rhombic sulfur S8 (s)';
                   FeS2_sediment_zt,    'Pyrite FeS2 (s)';
                   PO4_sediment_zt,     'Phosphate PO4(3-) (aq)';
                   PO4adsa_sediment_zt, 'Solid phosphorus pool a PO4adsa (s)';
                   PO4adsb_sediment_zt, 'Solid phosphorus pool b PO4adsb (s)';
                   NO3_sediment_zt,     'Nitrate NO3(-) (aq)';
                   NH4_sediment_zt,     'Ammonium NH4(+) (aq)';
                   Ca2_sediment_zt,     'Calcium Ca(2+) (aq)';
                   Ca3PO42_sediment_zt, 'Apatite Ca3PO42 (s)';
                   H_sediment_zt,       'H(+)(aq)';
                   OH_sediment_zt,      'OH(-)(aq)';
                   CO2_sediment_zt,     'CO2(aq)';
                   CO3_sediment_zt,     'CO3(2-)(aq)';
                   HCO3_sediment_zt,    'HCO3(-)(aq)';
                   NH3_sediment_zt,     'NH3(aq)';
                   H2CO3_sediment_zt,   'H2CO3(aq)';
                   pH_sediment_zt,      'pH in sediment';
                   OM_flux_sediment_zt, 'OM flux to sediment';
                   OM2_flux_sediment_zt,'OM2 flux to sediment';
                   O2_flux_sediment_zt, 'Oxygen flux WC to Sediments';
                   PO4_flux_sediment_zt,'PO4 flux WC to Sediments';
                   NO3_flux_sediment_zt,'NO3 flux to sediment';
                   FeOH3_flux_sediment_zt,'FeOH3 flux to sediment';
                   z_sediment,         'z';
                   % R_values_sediment_zt,'R values of TEA oxidations sediment';
                   Bioirrigation_sediment_zt, 'Fluxes of bioirrigation';
                   params_sediment,     'Sediments params';
};

x = linspace(0,100,128);
w = 1;
D = 368.53;
t = 1;
C_sol = 1/2*( erfc( (x-w*t)./2./sqrt(D*t) ) + exp(w.*x/D).*erfc((x+w.*t)/2/sqrt(D*t)));
plot(x, O2_sediment_zt(:,end), 'kx','MarkerSize', 12); hold on;plot(x, C_sol, 'LineWidth', 3);legend('numerical','analytical');

end

