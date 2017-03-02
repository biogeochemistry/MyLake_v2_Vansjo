function [ O2z, Pz, Fe2z, NO3z, NH4z ] = update_wc( O2z, Pz, Fe2z, NO3z, NH4z, MyLake_params, sediment_SWI_fluxes, sediment_bioirrigation_fluxes )
    %UPDATE_WC Summary of this function goes here

    global debugg

    % Diffusion Fluxes:
    O2z = update_C_due_to_flux(O2z, sediment_SWI_fluxes{1}, MyLake_params, true);
    Pz = update_C_due_to_flux(Pz, sediment_SWI_fluxes{4}, MyLake_params, true);
    Fe2z = update_C_due_to_flux(Fe2z, sediment_SWI_fluxes{7}, MyLake_params, true);
    NO3z = update_C_due_to_flux(NO3z, sediment_SWI_fluxes{5}, MyLake_params, true);
    NH4z = update_C_due_to_flux(NH4z, sediment_SWI_fluxes{8}, MyLake_params, true);

    % Boudreau, B.P., 1999. Metals and models : Diagenetic modelling in freshwater lacustrine sediments *. , pp.227–251.

    % Bioirrigation
    O2z = update_C_due_to_flux(O2z, sediment_bioirrigation_fluxes{1}, MyLake_params, true);
    Pz = update_C_due_to_flux(Pz, sediment_bioirrigation_fluxes{2}, MyLake_params, true);
    Fe2z = update_C_due_to_flux(Fe2z, sediment_bioirrigation_fluxes{3}, MyLake_params, true);
    NO3z = update_C_due_to_flux(NO3z, sediment_bioirrigation_fluxes{4}, MyLake_params, true);
    NH4z = update_C_due_to_flux(NH4z, sediment_bioirrigation_fluxes{5}, MyLake_params, true);


    if debugg
        if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
            error('NaN')
        end
    end
end

%% update_C_due_to_flux: Update concentration of the WC due to flux from/to sediment
function [C] = update_C_due_to_flux(C, flux, MyLake_params, is_sediment_source)
    Az_end = MyLake_params('Az(end)');
    Vz_end = MyLake_params('Vz(end)');
    dt     = MyLake_params('dt');

    if is_sediment_source
        dC = flux * dt * Az_end / Vz_end ;
    else
        dC = flux * dt * Az_end / Vz_end * (flux > 0);
    end
    C(end) = C(end) - dC;
    C(end) = (C(end) > 0).*C(end);
end


%
