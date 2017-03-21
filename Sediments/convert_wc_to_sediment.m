function [sediment_bc] = convert_wc_to_sediment(MyLake_concentrations, MyLake_params, sediment_params)
%convert_wc_to_sediment function convert BC values for sediment module (different Units)
    %   [Water-column units] ------>  [Sediments units]
    % for dissolved species - concentration BC
    % for solid - Neumann (flux) BC
    pH = MyLake_params.pH;
    % pH = 6.47;
    w_s = MyLake_params.w_s * 100 * 365; %settling velocity for S [m d-1] -> [cm year-1]
    w_chl = MyLake_params.w_chl * 100 * 365; %settling velocity for S [m d-1] -> [cm year-1]
    w_chl_2 = MyLake_params.w_chl_2 * 100 * 365; %settling velocity for S [m d-1] -> [cm year-1]
    fi = sediment_params.fi;

    sediment_bc.Ox_c = dissolved_bc(MyLake_concentrations.O2z);
    sediment_bc.OM1_fx = solid_bc(MyLake_concentrations.Chlz, w_chl, fi) + solid_bc(MyLake_concentrations.Cz, w_chl_2, fi);
    sediment_bc.OM2_fx = 0;
    sediment_bc.PO4_c = dissolved_bc(MyLake_concentrations.Pz);
    sediment_bc.NO3_c = dissolved_bc(MyLake_concentrations.NO3z);
    sediment_bc.FeOH3_fx = solid_bc(MyLake_concentrations.Fe3z, w_s, fi);
    sediment_bc.SO4_c = dissolved_bc(MyLake_concentrations.SO4z);
    sediment_bc.Fe2_c = dissolved_bc(MyLake_concentrations.Fe2z);
    sediment_bc.FeOOH_fx = 0;
    sediment_bc.FeS_fx = 0;
    sediment_bc.S0_c = 0;
    sediment_bc.S8_fx = 0;
    sediment_bc.FeS2_fx = 0;
    sediment_bc.AlOH3_fx = solid_bc(MyLake_concentrations.Al3z, w_s, fi);
    sediment_bc.PO4adsa_fx = solid_bc(MyLake_concentrations.PPz, w_s, fi);
    sediment_bc.PO4adsb_fx = 0;
    sediment_bc.Ca2_c = dissolved_bc(MyLake_concentrations.Ca2z);
    sediment_bc.Ca3PO42_fx = 0;
    sediment_bc.OMS_fx = 0;
    sediment_bc.H_c = 10^-pH*10^3;
    sediment_bc.OH_c = 10^-(14-pH)*10^3;
    sediment_bc.CO2_c = 0;
    sediment_bc.CO3_c = 2.19E-05;
    sediment_bc.HCO3_c = 0.62387047;
    sediment_bc.NH3_c = 3.68E-09;
    sediment_bc.NH4_c = dissolved_bc(MyLake_concentrations.NH4z);
    sediment_bc.HS_c = 1.01E-10;
    sediment_bc.H2S_c = 1.06E-10;
    sediment_bc.H2CO3_c = 1.06E-15;
    sediment_bc.DOM1_c = dissolved_bc(MyLake_concentrations.DOPz);
    sediment_bc.DOM2_c = dissolved_bc(MyLake_concentrations.DOCz);
    sediment_bc.T = MyLake_params.Tz(end);

if any(isnan(dissolved_bc(MyLake_concentrations.O2z))) | any(isnan(dissolved_bc(MyLake_concentrations.Pz))) | any(isnan(dissolved_bc(MyLake_concentrations.NO3z))) | any(isnan(dissolved_bc(MyLake_concentrations.SO4z))) | any(isnan(dissolved_bc(MyLake_concentrations.Fe2z))) | any(isnan(dissolved_bc(MyLake_concentrations.Ca2z))) | any(isnan(dissolved_bc(MyLake_concentrations.NH4z)))
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



