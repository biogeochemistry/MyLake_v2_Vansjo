function [sediment_bc] = convert_wc_to_sediment(MyLake_concentrations, MyLake_params, sediment_params)
%convert_wc_to_sediment function convert BC values for sediment module (different Units)
    %   [Water-column units] ------>  [Sediments units]
    % for dissolved species - concentration BC
    % for solid - Neumann (flux) BC
    pH = MyLake_params('pH');
    % pH = 6.47;
    w_s = MyLake_params('w_s') * 100 * 365; %settling velocity for Chl [m d-1] -> [cm year-1]
    fi = sediment_params('fi');


    sediment_bc = {...
        dissolved_bc(MyLake_concentrations('O2z')),                           'Ox_c';
        OM1_flux(MyLake_concentrations, MyLake_params, sediment_params),      'OM1_fx';
        % OM2_flux(MyLake_concentrations, MyLake_params),                       'OM2_fx';
        solid_bc(MyLake_concentrations('DOPz'), w_s/20, fi) + solid_bc(MyLake_concentrations('DOCz'), w_s/20, fi),                 'OM2_fx'; % NOTE:  We need to introduce DOC:POC ratio 20
        dissolved_bc(MyLake_concentrations('Pz')),                            'PO4_c';
        dissolved_bc(MyLake_concentrations('NO3z')),                          'NO3_c';
        solid_bc(MyLake_concentrations('Fe3z'), w_s, fi),                     'FeOH3_fx';
        dissolved_bc(MyLake_concentrations('SO4z')),                          'SO4_c';
        dissolved_bc(MyLake_concentrations('Fe2z')),                          'Fe2_c';
        0,                'FeOOH_fx'; % 0; %from Canavan et al AML
        0,                'FeS_fx'; % 0; % Flux for solid species
        0,                'S0_c'; % 0 ; % Exact concentration for solute species
        0,                'S8_fx'; % 0; %from Canavan et al AML
        0,                'FeS2_fx'; % 0; % Flux for solid species
        solid_bc(MyLake_concentrations('Al3z'), w_s, fi),                     'AlOH3_fx';
        solid_bc(MyLake_concentrations('PPz'), w_s, fi),                      'PO4adsa_fx';
        0,                'PO4adsb_fx'; % 0; % Flux for solid species
        dissolved_bc(MyLake_concentrations('Ca2z')),                          'Ca2_c';
        0,                'Ca3PO42_fx'; % 0; % Flux for solid species
        0,                'OMS_fx'; % 0; % Flux for solid species
        10^-pH*10^3,      'H_c';
        0.00007354365,    'OH_c';
        9.960074871,      'CO2_c';
        2.19E-05,         'CO3_c';
        0.62387047,       'HCO3_c';
        3.68E-09,         'NH3_c';
        dissolved_bc(MyLake_concentrations('NH4z')),                           'NH4_c';
        1.01E-10,         'HS_c';
        1.06E-10,         'H2S_c';
        1,      'H2CO3_c';
    };

    % sediment_bc = {...
    %         0  'Ox_c';
    %         0, 'OM1_fx';
    %         0, 'OM2_fx';
    %         0, 'PO4_c';
    %         0, 'NO3_c';
    %         0, 'FeOH3_fx';
    %         0, 'SO4_c';
    %         0, 'Fe2_c';
    %         0, 'FeOOH_fx'; % 0; %from Canavan et al AML
    %         0, 'FeS_fx'; % 0; % Flux for solid species
    %         0, 'S0_c'; % 0 ; % Exact concentration for solute species
    %         0, 'S8_fx'; % 0; %from Canavan et al AML
    %         0, 'FeS2_fx'; % 0; % Flux for solid species
    %         0, 'AlOH3_fx';
    %         0, 'PO4adsa_fx';
    %         0, 'PO4adsb_fx'; % 0; % Flux for solid species
    %         0, 'Ca2_c';
    %         0, 'Ca3PO42_fx'; % 0; % Flux for solid species
    %         0, 'OMS_fx'; % 0; % Flux for solid species
    %         0, 'H_c';
    %         0, 'OH_c';
    %         0, 'CO2_c';
    %         0, 'CO3_c';
    %         0, 'HCO3_c';
    %         0, 'NH3_c';
    %         0, 'NH4_c';
    %         0, 'HS_c';
    %         0, 'H2S_c';
    %         0, 'H2CO3_c';
    % };
    sediment_bc = containers.Map({sediment_bc{:,2}},{sediment_bc{:,1}});

if any(isnan(dissolved_bc(MyLake_concentrations('O2z')))) | any(isnan(dissolved_bc(MyLake_concentrations('Pz')))) | any(isnan(dissolved_bc(MyLake_concentrations('NO3z')))) | any(isnan(dissolved_bc(MyLake_concentrations('SO4z')))) | any(isnan(dissolved_bc(MyLake_concentrations('Fe2z')))) | any(isnan(dissolved_bc(MyLake_concentrations('Ca2z')))) | any(isnan(dissolved_bc(MyLake_concentrations('NH4z'))))
    error('stop')
end

end

function C_bc = dissolved_bc(C)
% return the value of boundary concentration for sediment
% C - concentration of the particular species in MyLake [umol/cm3]
    C_bc = C(end);
end

function solid_fx = solid_bc(C, w_s, fi)
    % C   - concentration in WC [mg m-3]
    % M_C - molar mass    [mg mol-1]
    % w_s - settling velocity of solids [cm year-1]
    % fi  - porosity [-]
    % solid_fx - flux of solid at SWI [umol cm-2 yr-1]
    solid_fx = (1 - fi(1)) * w_s  * C(end);
end

function OM2_fx = OM2_flux(MyLake_concentrations, MyLake_params)
    % this function convert H_netsed_catch to OM2 fluxes
    % H_netsed_catch - flux of carbon from catchments (profile) [m day-1, dry]
    % SS_C - Carbon fraction in H_netsed_catch [~]
    % density_org_H_nc - Density of organic fraction in H_netsed_catch [g cm-3]
    % OM2_fx - Sediments flux of OM2 [umol cm-2 yr-1]
    H_netsed_catch = MyLake_concentrations('H_netsed_catch');
    SS_C = MyLake_params('SS_C');
    density_org_H_nc = MyLake_params('density_org_H_nc');
    M_C = 12*1e-6; % Molar weight of carbon (g umol-1)
    OM2_fx = H_netsed_catch(end)*100*365*SS_C*density_org_H_nc/M_C; %
end

function OM1_fx = OM1_flux(MyLake_concentrations, MyLake_params, sediment_params)
% function convert flux of Chl to OM1 flux in sediment
% Chl - chlorophyll (group 1) [mg m-3]
% Cz - chlorophyll (group 2) [mg m-3]
% w_chl - settling velocity for Chl a (m day-1)
% Cx1 -
% fi -
% M_C - Molar mass of carbon [mg mol-1]
% Mass_Ratio_C_Chl - Fixed empirical ratio C:Chl (mass/mass)

    Chl = MyLake_concentrations('Chlz');
    Cz = MyLake_concentrations('Cz');
    w_chl = MyLake_params('w_chl') * 100 * 365; %settling velocity for Chl [m d-1] -> [cm year-1]
    Cx1= sediment_params('Cx1');
    fi = sediment_params('fi');
    M_C = 12011; % Molar mass of carbon [mg mol-1]
    Mass_Ratio_C_Chl = MyLake_params('Mass_Ratio_C_Chl');

    OM1_fx =  (1 - fi(1)) * w_chl * (Chl(end)+Cz(end)) * Mass_Ratio_C_Chl / ( M_C * Cx1);
end


