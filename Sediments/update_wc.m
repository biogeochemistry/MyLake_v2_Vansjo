function [ O2z, Pz, Fe2z, NO3z, NH4z ] = update_wc( O2z, Pz, Fe2z, NO3z, NH4z, MyLake_params, sediment_SWI_fluxes, sediment_bioirrigation_fluxes )
    %UPDATE_WC Summary of this function goes here

    % Diffusion Fluxes:
    O2z = update_C_due_to_flux(O2z, sediment_SWI_fluxes.Ox, MyLake_params);
    Pz = update_C_due_to_flux(Pz, sediment_SWI_fluxes.PO4, MyLake_params);
    Fe2z = update_C_due_to_flux(Fe2z, sediment_SWI_fluxes.Fe2, MyLake_params);
    NO3z = update_C_due_to_flux(NO3z, sediment_SWI_fluxes.NO3, MyLake_params);
    NH4z = update_C_due_to_flux(NH4z, sediment_SWI_fluxes.NH4, MyLake_params);

    % % Boudreau, B.P., 1999. Metals and models : Diagenetic modelling in freshwater lacustrine sediments *. , pp.227–251.

    % % Bioirrigation
    O2z = update_C_due_to_flux(O2z, sediment_bioirrigation_fluxes.Ox, MyLake_params);
    Pz = update_C_due_to_flux(Pz, sediment_bioirrigation_fluxes.PO4, MyLake_params);
    Fe2z = update_C_due_to_flux(Fe2z, sediment_bioirrigation_fluxes.Fe2, MyLake_params);
    NO3z = update_C_due_to_flux(NO3z, sediment_bioirrigation_fluxes.NO3, MyLake_params);
    NH4z = update_C_due_to_flux(NH4z, sediment_bioirrigation_fluxes.NH4, MyLake_params);



    if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
        error('NaN')
    end
end

%% update_C_due_to_flux: Update concentration of the WC due to flux from/to sediment
function [C] = update_C_due_to_flux(C, flux, MyLake_params)
    Az_end = MyLake_params.Az(end);
    Vz_end = MyLake_params.Vz(end);
    dt     = MyLake_params.dt;

    dC = flux * dt * Az_end / Vz_end;
    C(end) = C(end) + dC;
    C(end) = (C(end) > 0).*C(end);
end


%
