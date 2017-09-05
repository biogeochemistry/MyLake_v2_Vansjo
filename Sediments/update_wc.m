function [mylake_new_resutls] = update_wc( mylake_prev_results, sediment_concentrations, sediment_transport_fluxes, sediment_bioirrigation_fluxes, MyLake_params, sediment_params, effective_depth)
    %UPDATE_WC Summary of this function goes here

    effective_depth(isnan(effective_depth)) = 0;

    O2z = mylake_prev_results.O2z;
    Pz = mylake_prev_results.Pz;
    Fe2z = mylake_prev_results.Fe2z;
    NO3z = mylake_prev_results.NO3z;
    NH4z = mylake_prev_results.NH4z;
    SO4z = mylake_prev_results.SO4z;
    HSz = mylake_prev_results.HSz;
    DOCz = mylake_prev_results.DOCz;
    DOPz = mylake_prev_results.DOPz;
    CH4aqz = mylake_prev_results.CH4aqz;
    CH4gz = mylake_prev_results.CH4gz;
    CO2z = mylake_prev_results.CO2z;
    HCO3z = mylake_prev_results.HCO3z;
    CO3z = mylake_prev_results.CO3z;
    Chlz= mylake_prev_results.Chlz;
    Cz= mylake_prev_results.Cz;
    POPz= mylake_prev_results.POPz;
    POCz= mylake_prev_results.POCz;
    Fe3z= mylake_prev_results.Fe3z;
    FeSz= mylake_prev_results.FeSz;
    Al3z= mylake_prev_results.Al3z;
    PPz= mylake_prev_results.PPz;
    CaCO3z= mylake_prev_results.CaCO3z;

    % Transport Fluxes dissolved:
    O2z = update_bottom_concentrations_below_eff_depth_aq(O2z, sediment_transport_fluxes.O2, MyLake_params, sediment_params, effective_depth);
    Pz = update_bottom_concentrations_below_eff_depth_aq(Pz, sediment_transport_fluxes.PO4, MyLake_params, sediment_params, effective_depth);
    Fe2z = update_bottom_concentrations_below_eff_depth_aq(Fe2z, sediment_transport_fluxes.Fe2, MyLake_params, sediment_params, effective_depth);
    NO3z = update_bottom_concentrations_below_eff_depth_aq(NO3z, sediment_transport_fluxes.NO3, MyLake_params, sediment_params, effective_depth);
    NH4z = update_bottom_concentrations_below_eff_depth_aq(NH4z, sediment_transport_fluxes.NH4, MyLake_params, sediment_params, effective_depth);
    SO4z = update_bottom_concentrations_below_eff_depth_aq(SO4z, sediment_transport_fluxes.SO4, MyLake_params, sediment_params, effective_depth);
    HSz = update_bottom_concentrations_below_eff_depth_aq(HSz, sediment_transport_fluxes.HS, MyLake_params, sediment_params, effective_depth);
    DOPz = update_bottom_concentrations_below_eff_depth_aq(DOPz, sediment_transport_fluxes.DOP, MyLake_params, sediment_params, effective_depth);
    DOCz = update_bottom_concentrations_below_eff_depth_aq(DOCz, sediment_transport_fluxes.DOC, MyLake_params, sediment_params, effective_depth);
    CH4aqz = update_bottom_concentrations_below_eff_depth_aq(CH4aqz, sediment_transport_fluxes.CH4aq, MyLake_params, sediment_params, effective_depth);
    CH4gz = update_bottom_concentrations_below_eff_depth_aq(CH4gz, sediment_transport_fluxes.CH4g, MyLake_params, sediment_params, effective_depth);
    CO2z = update_bottom_concentrations_below_eff_depth_aq(CO2z, sediment_transport_fluxes.CO2, MyLake_params, sediment_params, effective_depth);
    HCO3z = update_bottom_concentrations_below_eff_depth_aq(HCO3z, sediment_transport_fluxes.HCO3, MyLake_params, sediment_params, effective_depth);
    CO3z = update_bottom_concentrations_below_eff_depth_aq(CO3z, sediment_transport_fluxes.CO3, MyLake_params, sediment_params, effective_depth);

    % Transport Fluxes solid:
    % Chlz = update_bottom_concentrations_below_eff_depth_solid(Chlz, Chlz./(Chlz+Cz+1e-16)*sediment_transport_fluxes.Chl, MyLake_params, sediment_params, effective_depth);
    % Cz = update_bottom_concentrations_below_eff_depth_solid(Cz, Cz./(Chlz+Cz+1e-16)*sediment_transport_fluxes.Chl, MyLake_params, sediment_params, effective_depth);
    % POPz = update_bottom_concentrations_below_eff_depth_solid(POPz, sediment_transport_fluxes.POP, MyLake_params, sediment_params, effective_depth);
    % POCz = update_bottom_concentrations_below_eff_depth_solid(POCz, sediment_transport_fluxes.POC, MyLake_params, sediment_params, effective_depth);
    % Fe3z = update_bottom_concentrations_below_eff_depth_solid(Fe3z, sediment_transport_fluxes.FeOH3, MyLake_params, sediment_params, effective_depth);
    % FeSz = update_bottom_concentrations_below_eff_depth_solid(FeSz, sediment_transport_fluxes.FeS, MyLake_params, sediment_params, effective_depth);
    % Al3z = update_bottom_concentrations_below_eff_depth_solid(Al3z, sediment_transport_fluxes.AlOH3, MyLake_params, sediment_params, effective_depth);
    % PPz = update_bottom_concentrations_below_eff_depth_solid(PPz, sediment_transport_fluxes.PO4adsa, MyLake_params, sediment_params, effective_depth);
    % CaCO3z = update_bottom_concentrations_below_eff_depth_solid(CaCO3z, sediment_transport_fluxes.CaCO3, MyLake_params, sediment_params, effective_depth);


    % % Boudreau, B.P., 1999. Metals and models : Diagenetic modelling in freshwater lacustrine sediments *. , pp.227–251.

    % % Bioirrigation dissolved only
    O2z = update_bottom_concentrations_below_eff_depth_aq(O2z, sediment_bioirrigation_fluxes.O2, MyLake_params, sediment_params, effective_depth);
    Pz = update_bottom_concentrations_below_eff_depth_aq(Pz, sediment_bioirrigation_fluxes.PO4, MyLake_params, sediment_params, effective_depth);
    Fe2z = update_bottom_concentrations_below_eff_depth_aq(Fe2z, sediment_bioirrigation_fluxes.Fe2, MyLake_params, sediment_params, effective_depth);
    NO3z = update_bottom_concentrations_below_eff_depth_aq(NO3z, sediment_bioirrigation_fluxes.NO3, MyLake_params, sediment_params, effective_depth);
    NH4z = update_bottom_concentrations_below_eff_depth_aq(NH4z, sediment_bioirrigation_fluxes.NH4, MyLake_params, sediment_params, effective_depth);
    SO4z = update_bottom_concentrations_below_eff_depth_aq(SO4z, sediment_bioirrigation_fluxes.SO4, MyLake_params, sediment_params, effective_depth);
    HSz = update_bottom_concentrations_below_eff_depth_aq(HSz, sediment_bioirrigation_fluxes.HS, MyLake_params, sediment_params, effective_depth);
    DOPz = update_bottom_concentrations_below_eff_depth_aq(DOPz, sediment_bioirrigation_fluxes.DOP, MyLake_params, sediment_params, effective_depth);
    DOCz = update_bottom_concentrations_below_eff_depth_aq(DOCz, sediment_bioirrigation_fluxes.DOC, MyLake_params, sediment_params, effective_depth);
    CH4aqz = update_bottom_concentrations_below_eff_depth_aq(CH4aqz, sediment_bioirrigation_fluxes.CH4aq, MyLake_params, sediment_params, effective_depth);
    CH4gz = update_bottom_concentrations_below_eff_depth_aq(CH4gz, sediment_bioirrigation_fluxes.CH4g, MyLake_params, sediment_params, effective_depth);
    CO2z = update_bottom_concentrations_below_eff_depth_aq(CO2z, sediment_bioirrigation_fluxes.CO2, MyLake_params, sediment_params, effective_depth);
    HCO3z = update_bottom_concentrations_below_eff_depth_aq(HCO3z, sediment_bioirrigation_fluxes.HCO3, MyLake_params, sediment_params, effective_depth);
    CO3z = update_bottom_concentrations_below_eff_depth_aq(CO3z, sediment_bioirrigation_fluxes.CO3, MyLake_params, sediment_params, effective_depth);

    mylake_new_resutls.O2z = O2z;
    mylake_new_resutls.Pz = Pz;
    mylake_new_resutls.Fe2z = Fe2z;
    mylake_new_resutls.NO3z = NO3z;
    mylake_new_resutls.NH4z = NH4z;
    mylake_new_resutls.SO4z = SO4z;
    mylake_new_resutls.HSz = HSz;
    mylake_new_resutls.DOPz = DOPz;
    mylake_new_resutls.DOCz = DOCz;
    mylake_new_resutls.CH4aqz = CH4aqz;
    mylake_new_resutls.CH4gz = CH4gz;
    mylake_new_resutls.CO2z = CO2z;
    mylake_new_resutls.HCO3z = HCO3z;
    mylake_new_resutls.CO3z = CO3z;
    mylake_new_resutls.Chlz = Chlz;
    mylake_new_resutls.Cz = Cz;
    mylake_new_resutls.POPz = POPz;
    mylake_new_resutls.POCz = POCz;
    mylake_new_resutls.Fe3z = Fe3z;
    mylake_new_resutls.FeSz = FeSz;
    mylake_new_resutls.Al3z = Al3z;
    mylake_new_resutls.PPz = PPz;
    mylake_new_resutls.CaCO3z = CaCO3z;

    if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
        error('NaN')
    end
end


%% update_bottom_concentrations_below_eff_depth: Update concentration of the WC due to flux from/to sediment bellow effective depth
function [C] = update_bottom_concentrations_below_eff_depth_aq(C, flux, MyLake_params, sediment_params, effective_depth)

    Az_dif = -diff([MyLake_params.Az; 0]);
    Vz = MyLake_params.Vz;
    dt     = MyLake_params.dt;
    % effective_depth = sediment_params.effective_depth;
    z = MyLake_params.zz;

    dC = flux * dt * Az_dif ./ Vz;
    C = C + (z > effective_depth).*dC;
    C = (C > 0).*C;
end

%% update_bottom_concentrations_below_eff_depth: Update concentration of the WC due to flux to sediment bellow effective depth
function [C] = update_bottom_concentrations_below_eff_depth_solid(C, flux, MyLake_params, sediment_params, effective_depth)

    Az_dif = -diff([MyLake_params.Az; 0]);
    Az_dif(end) = 0;
    Vz = MyLake_params.Vz;
    dt     = MyLake_params.dt;
    % effective_depth = sediment_params.effective_depth;
    z = MyLake_params.zz;

    dC = flux .* dt .* Az_dif ./ Vz;
    C = C + (z > effective_depth).*dC;
    C = (C > 0).*C;
end


%% update_only_bottom_concentrations: Update concentration of the WC due to flux from/to sediment
function [C] = update_only_bottom_concentrations(C, flux, MyLake_params, sediment_params)
    Az_end = MyLake_params.Az(end);
    Vz_end = MyLake_params.Vz(end);
    dt     = MyLake_params.dt;
    phi = sediment_params.phi;


    dC = flux * dt * Az_end / Vz_end;
    C(end) = C(end) + dC;
    C(end) = (C(end) > 0).*C(end);
end

%% update_C_as_dirichlet: function description
function [C_wc] = update_C_as_dirichlet(C_wc, BC_value)
    % BC_value - boundary condition produced by sediment
    % C_wc - concentration of the species in the water-column
    C_wc(end) = BC_value;
end
