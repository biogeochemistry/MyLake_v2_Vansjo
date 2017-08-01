function [mylake_new_resutls] = update_wc( mylake_prev_results, sediment_concentrations, sediment_transport_fluxes, sediment_bioirrigation_fluxes, MyLake_params, sediment_params)
    %UPDATE_WC Summary of this function goes here

    O2z = mylake_prev_results.O2z;
    Pz = mylake_prev_results.Pz;
    Fe2z = mylake_prev_results.Fe2z;
    NO3z = mylake_prev_results.NO3z;
    NH4z = mylake_prev_results.NH4z;
    DOCz = mylake_prev_results.DOCz;
    DOPz = mylake_prev_results.DOPz;

    % Diffusion Fluxes:
    O2z = update_C_as_neumann(O2z, sediment_transport_fluxes.O2, MyLake_params, sediment_params);
    Pz = update_C_as_neumann(Pz, sediment_transport_fluxes.PO4, MyLake_params, sediment_params);
    Fe2z = update_C_as_neumann(Fe2z, sediment_transport_fluxes.Fe2, MyLake_params, sediment_params);
    NO3z = update_C_as_neumann(NO3z, sediment_transport_fluxes.NO3, MyLake_params, sediment_params);
    NH4z = update_C_as_neumann(NH4z, sediment_transport_fluxes.NH4, MyLake_params, sediment_params);
    DOPz = update_C_as_neumann(DOPz, sediment_transport_fluxes.DOP, MyLake_params, sediment_params);
    DOCz = update_C_as_neumann(DOCz, sediment_transport_fluxes.DOC, MyLake_params, sediment_params);

    % % Boudreau, B.P., 1999. Metals and models : Diagenetic modelling in freshwater lacustrine sediments *. , pp.227–251.

    % % Bioirrigation
    O2z = update_C_as_neumann(O2z, sediment_bioirrigation_fluxes.O2, MyLake_params, sediment_params);
    Pz = update_C_as_neumann(Pz, sediment_bioirrigation_fluxes.PO4, MyLake_params, sediment_params);
    Fe2z = update_C_as_neumann(Fe2z, sediment_bioirrigation_fluxes.Fe2, MyLake_params, sediment_params);
    NO3z = update_C_as_neumann(NO3z, sediment_bioirrigation_fluxes.NO3, MyLake_params, sediment_params);
    NH4z = update_C_as_neumann(NH4z, sediment_bioirrigation_fluxes.NH4, MyLake_params, sediment_params);
    DOPz = update_C_as_neumann(DOPz, sediment_bioirrigation_fluxes.DOP, MyLake_params, sediment_params);
    DOCz = update_C_as_neumann(DOCz, sediment_bioirrigation_fluxes.DOC, MyLake_params, sediment_params);

    mylake_new_resutls.O2z = O2z;
    mylake_new_resutls.Pz = Pz;
    mylake_new_resutls.Fe2z = Fe2z;
    mylake_new_resutls.NO3z = NO3z;
    mylake_new_resutls.NH4z = NH4z;
    mylake_new_resutls.DOPz = DOPz;
    mylake_new_resutls.DOCz = DOCz;

    if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
        error('NaN')
    end
end

%% update_C_as_neumann: Update concentration of the WC due to flux from/to sediment
function [C] = update_C_as_neumann(C, flux, MyLake_params, sediment_params)
    Az_end = MyLake_params.Az(end);
    Vz_end = MyLake_params.Vz(end);
    dt     = MyLake_params.dt;
    phi = sediment_params.phi;

    dC = flux * dt * Az_end / Vz_end;
    C(end) = C(end) + dC;
    C(end) = (C(end) > 0).*C(end);
end

%% update_C_as_dirichlet: function description
function [C] = update_C_as_dirichlet(C_wc, BC_value)
    % BC_value - boundary condition produced by sediment
    % C_wc - concentration of the species in the water-column
    C_WC(end) = BC_value;
end
