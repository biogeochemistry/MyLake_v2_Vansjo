function [sediment_bc] = update_sediment(mylake_temp_results, mylake_params, sediment_params, effective_depth)
%update_sedimets function convert BC values for sediment module (different Units)
    %   [Water-column units] ------>  [Sediments units]
    % for dissolved species - concentration BC
    % for solid - Neumann (flux) BC

    % pH = 6.47;
    w_s = mylake_params.w_s * 100 * 365; %settling velocity for S [m d-1] -> [cm year-1]
    w_chl = mylake_params.w_chl * 100 * 365; %settling velocity for S [m d-1] -> [cm year-1]
    w_chl_2 = mylake_params.w_chl_2 * 100 * 365; %settling velocity for S [m d-1] -> [cm year-1]
    phi = sediment_params.phi;

    effective_depth(isnan(effective_depth)) = 0;

    sediment_bc.O2_c = dissolved_bc_with_effective_depth(mylake_temp_results.O2z, sediment_params, mylake_params, effective_depth);
    sediment_bc.Chl_fx = solid_bc_with_effective_depth(mylake_temp_results.Chlz, w_chl, sediment_params, mylake_params, effective_depth) + solid_bc_with_effective_depth(mylake_temp_results.Cz, w_chl_2, sediment_params, mylake_params, effective_depth);
    sediment_bc.POP_fx = solid_bc_with_effective_depth(mylake_temp_results.POPz, w_s, sediment_params, mylake_params, effective_depth);
    sediment_bc.POC_fx = solid_bc_with_effective_depth(mylake_temp_results.POCz, w_s, sediment_params, mylake_params, effective_depth);
    sediment_bc.PO4_c = dissolved_bc_with_effective_depth(mylake_temp_results.Pz, sediment_params, mylake_params, effective_depth);
    sediment_bc.NO3_c = dissolved_bc_with_effective_depth(mylake_temp_results.NO3z, sediment_params, mylake_params, effective_depth);
    sediment_bc.FeOH3_fx = solid_bc_with_effective_depth(mylake_temp_results.Fe3z, w_s, sediment_params, mylake_params, effective_depth);
    sediment_bc.SO4_c = dissolved_bc_with_effective_depth(mylake_temp_results.SO4z, sediment_params, mylake_params, effective_depth);
    sediment_bc.Fe2_c = dissolved_bc_with_effective_depth(mylake_temp_results.Fe2z, sediment_params, mylake_params, effective_depth);
    sediment_bc.FeOOH_fx = 0;
    sediment_bc.FeS_fx = solid_bc_with_effective_depth(mylake_temp_results.FeSz, w_s, sediment_params, mylake_params, effective_depth);
    sediment_bc.S0_c = 0;
    sediment_bc.S8_fx = 0;
    sediment_bc.FeS2_fx = 0;
    sediment_bc.AlOH3_fx = solid_bc_with_effective_depth(mylake_temp_results.Al3z, w_s, sediment_params, mylake_params, effective_depth);
    sediment_bc.PO4adsa_fx = solid_bc_with_effective_depth(mylake_temp_results.PPz, w_s, sediment_params, mylake_params, effective_depth);
    sediment_bc.PO4adsb_fx = 0;
    sediment_bc.Ca2_c = dissolved_bc_with_effective_depth(mylake_temp_results.Ca2z, sediment_params, mylake_params, effective_depth);
    sediment_bc.Ca3PO42_fx = 0;
    sediment_bc.OMS_fx = 0;
    sediment_bc.H3O_c = 10.^-7.5*10^3; % 10.^-mylake_temp_results.pHz(end)*10^3;
    sediment_bc.CaCO3_fx = solid_bc_with_effective_depth(mylake_temp_results.CaCO3z, w_s, sediment_params, mylake_params, effective_depth);
    sediment_bc.CO2_c = dissolved_bc_with_effective_depth(mylake_temp_results.CO2z, sediment_params, mylake_params, effective_depth); % gas
    sediment_bc.CO3_c = dissolved_bc_with_effective_depth(mylake_temp_results.CO3z, sediment_params, mylake_params, effective_depth);;
    sediment_bc.HCO3_c = dissolved_bc_with_effective_depth(mylake_temp_results.HCO3z, sediment_params, mylake_params, effective_depth);;
    sediment_bc.CO2g_c = 0;
    sediment_bc.NH3_c = 0;
    sediment_bc.NH4_c = dissolved_bc_with_effective_depth(mylake_temp_results.NH4z, sediment_params, mylake_params, effective_depth);
    sediment_bc.HS_c = dissolved_bc_with_effective_depth(mylake_temp_results.HSz, sediment_params, mylake_params, effective_depth);
    sediment_bc.H2S_c = 0;
    sediment_bc.DOP_c = dissolved_bc_with_effective_depth(mylake_temp_results.DOPz, sediment_params, mylake_params, effective_depth);
    sediment_bc.DOC_c = dissolved_bc_with_effective_depth(mylake_temp_results.DOCz, sediment_params, mylake_params, effective_depth);
    sediment_bc.CH4aq_c = dissolved_bc_with_effective_depth(mylake_temp_results.CH4aqz, sediment_params, mylake_params, effective_depth);
    sediment_bc.CH4g_fx = 0; % dissolved_bc_with_effective_depth(mylake_temp_results.CH4gz, sediment_params, mylake_params, effective_depth);
    sediment_bc.FeCO3_fx = 0; % dissolved_bc_with_effective_depth(mylake_temp_results.CH4gz, sediment_params, mylake_params, effective_depth);
    sediment_bc.Fe3PO42_fx = 0; % dissolved_bc_with_effective_depth(mylake_temp_results.CH4gz, sediment_params, mylake_params, effective_depth);
    sediment_bc.T = mylake_temp_results.Tz(end);

end

function C_bc = dissolved_bc_with_effective_depth(C, sediment_params, mylake_params, effective_depth)

    z = mylake_params.zz;
    % effective_depth = sediment_params.effective_depth;
    Az_dif = -diff([mylake_params.Az; 0]) .* (z > effective_depth);

    % estimate area average concentration at the bottom
    C_bc = sum(C .* Az_dif) / sum(Az_dif);
end



%% solid_bc_with_effective_depth: estimates average flux under effective sediment depth
function solid_fx = solid_bc_with_effective_depth(C, w_s, sediment_params, mylake_params, effective_depth)
    z = mylake_params.zz;
    % effective_depth = sediment_params.effective_depth;
    Az_dif = -diff([mylake_params.Az; 0]) .* (z > effective_depth);

    % estimate area average flux at the bottom
    solid_fx = w_s * sum(C .* Az_dif) / sum(Az_dif);
end




function C_bc = dissoved_bc(C, sediment_params, mylake_params)
% return the value of boundary concentration for sediment
% In MyLake the concentrations are bulk concentrations
% in sediment the are per V of H2O or per V of solid
% C - concentration of the particular species in MyLake [umol/cm3]
    C_bc = C(end);
end



function solid_fx = solid_bc(C, w_s, sediment_params, mylake_params)
    % C   - concentration in WC [umol/cm3]
    % w_s - settling velocity of solids [cm year-1]
    % mmol arriving at SWI interface is
    % solid_fx - flux of solid at SWI [umol cm-2 yr-1]

    solid_fx = w_s * C(end); % * (1 - phi(1));
end




