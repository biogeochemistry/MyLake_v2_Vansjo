function [MyLakeNewConcentrations] = update_wc( MyLakeOldConcentrations, MyLake_params, sediment_params,sediment_transport_fluxes, sediment_bioirrigation_fluxes )
    %UPDATE_WC Summary of this function goes here

    O2z = MyLakeOldConcentrations.O2z;
    Pz = MyLakeOldConcentrations.Pz;
    Fe2z = MyLakeOldConcentrations.Fe2z;
    NO3z = MyLakeOldConcentrations.NO3z;
    NH4z = MyLakeOldConcentrations.NH4z;
    DOCz = MyLakeOldConcentrations.DOCz;
    DOPz = MyLakeOldConcentrations.DOPz;

    % Diffusion Fluxes:
    O2z = update_C_due_to_flux(O2z, sediment_transport_fluxes.Ox, MyLake_params, sediment_params);
    Pz = update_C_due_to_flux(Pz, sediment_transport_fluxes.PO4, MyLake_params, sediment_params);
    Fe2z = update_C_due_to_flux(Fe2z, sediment_transport_fluxes.Fe2, MyLake_params, sediment_params);
    NO3z = update_C_due_to_flux(NO3z, sediment_transport_fluxes.NO3, MyLake_params, sediment_params);
    NH4z = update_C_due_to_flux(NH4z, sediment_transport_fluxes.NH4, MyLake_params, sediment_params);
    DOPz = update_C_due_to_flux(DOPz, sediment_transport_fluxes.DOM1, MyLake_params, sediment_params);
    DOCz = update_C_due_to_flux(DOCz, sediment_transport_fluxes.DOM2, MyLake_params, sediment_params);

    % % Boudreau, B.P., 1999. Metals and models : Diagenetic modelling in freshwater lacustrine sediments *. , pp.227–251.

    % % Bioirrigation
    O2z = update_C_due_to_flux(O2z, sediment_bioirrigation_fluxes.Ox, MyLake_params, sediment_params);
    Pz = update_C_due_to_flux(Pz, sediment_bioirrigation_fluxes.PO4, MyLake_params, sediment_params);
    Fe2z = update_C_due_to_flux(Fe2z, sediment_bioirrigation_fluxes.Fe2, MyLake_params, sediment_params);
    NO3z = update_C_due_to_flux(NO3z, sediment_bioirrigation_fluxes.NO3, MyLake_params, sediment_params);
    NH4z = update_C_due_to_flux(NH4z, sediment_bioirrigation_fluxes.NH4, MyLake_params, sediment_params);
    DOPz = update_C_due_to_flux(DOPz, sediment_bioirrigation_fluxes.DOM1, MyLake_params, sediment_params);
    DOCz = update_C_due_to_flux(DOCz, sediment_bioirrigation_fluxes.DOM2, MyLake_params, sediment_params);

    MyLakeNewConcentrations.O2z = O2z;
    MyLakeNewConcentrations.Pz = Pz;
    MyLakeNewConcentrations.Fe2z = Fe2z;
    MyLakeNewConcentrations.NO3z = NO3z;
    MyLakeNewConcentrations.NH4z = NH4z;
    MyLakeNewConcentrations.DOPz = DOPz;
    MyLakeNewConcentrations.DOCz = DOCz;

    if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
        error('NaN')
    end
end

%% update_C_due_to_flux: Update concentration of the WC due to flux from/to sediment
function [C] = update_C_due_to_flux(C, flux, MyLake_params, sediment_params)
    Az_end = MyLake_params.Az(end);
    Vz_end = MyLake_params.Vz(end);
    dt     = MyLake_params.dt;
    fi = sediment_params.fi;

    dC = fi(1)*flux * dt * Az_end / Vz_end;
    C(end) = C(end) + dC;
    C(end) = (C(end) > 0).*C(end);
end


%
