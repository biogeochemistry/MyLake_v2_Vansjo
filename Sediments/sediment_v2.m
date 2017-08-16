function [ sediment_bioirrigation_fluxes, sediment_transport_fluxes, sediment_concentrations, sediment_additional_results] = sediment_v2(sediment_concentrations, sediment_params, sediment_matrix_templates, sediment_bc)
  % SEDIMENTS This function models the chemical process in the sediment


  O2(:,1) = sediment_concentrations.O2;
  POP(:,1) = sediment_concentrations.POP;
  POC(:,1) = sediment_concentrations.POC;
  NO3(:,1) = sediment_concentrations.NO3;
  FeOH3(:,1) = sediment_concentrations.FeOH3;
  SO4(:,1) = sediment_concentrations.SO4;
  Fe2(:,1) = sediment_concentrations.Fe2;
  FeOOH(:,1) = sediment_concentrations.FeOOH;
  FeS(:,1) = sediment_concentrations.FeS;
  S0(:,1) = sediment_concentrations.S0;
  PO4(:,1)  = sediment_concentrations.PO4;
  S8(:,1) = sediment_concentrations.S8;
  FeS2(:,1)  = sediment_concentrations.FeS2;
  AlOH3(:,1) = sediment_concentrations.AlOH3;
  PO4adsa(:,1) = sediment_concentrations.PO4adsa;
  PO4adsb(:,1) = sediment_concentrations.PO4adsb;
  Ca2(:,1) = sediment_concentrations.Ca2;
  Ca3PO42(:,1) = sediment_concentrations.Ca3PO42;
  OMS(:,1) = sediment_concentrations.OMS;
  H(:,1) = sediment_concentrations.H;
  OH(:,1) = sediment_concentrations.OH;
  CO2(:,1) = sediment_concentrations.CO2;
  CO3(:,1) = sediment_concentrations.CO3;
  HCO3(:,1) = sediment_concentrations.HCO3;
  NH3(:,1) = sediment_concentrations.NH3;
  NH4(:,1) = sediment_concentrations.NH4;
  HS(:,1) = sediment_concentrations.HS;
  H2S(:,1) = sediment_concentrations.H2S;
  H2CO3(:,1) = sediment_concentrations.H2CO3;
  DOP(:,1) = sediment_concentrations.DOP;
  DOC(:,1) = sediment_concentrations.DOC;
  Chl(:,1) = sediment_concentrations.Chl;
  CH4(:,1) = sediment_concentrations.CH4;


  % model domain:
  years = sediment_params.years; %1 day
  ts    = sediment_params.ts; % time step
  alfax = sediment_params.alfax;
  D_O2  = sediment_params.D_O2;
  D_NO3 = sediment_params.D_NO3;
  D_SO4 = sediment_params.D_SO4;
  D_NH4 = sediment_params.D_NH4;
  D_Fe2 = sediment_params.D_Fe2;
  D_H2S = sediment_params.D_H2S;
  D_S0  = sediment_params.D_S0;
  D_PO4 = sediment_params.D_PO4;
  D_Ca2 = sediment_params.D_Ca2;
  D_HS  = sediment_params.D_HS;
  D_DOP = sediment_params.D_DOP;
  D_DOC = sediment_params.D_DOC;
  D_CH4 = sediment_params.D_DOC;
  phi    = sediment_params.phi;
  x  = sediment_params.x;

  dx = x(2)-x(1);
  xz = x'/100;
  % time domain:
  t     = 0:ts:years; % years
  m     = size(t,2); %steps in time
  dt    = t(2)-t(1);



  % Allocation of the memory and formation of template matrix:
  % =======================================================================================================

  % Solid species: row #1 of the cell "sediment_matrix_templates" is the solid template matrix
  [POP_AL, POP_AR]= sediment_matrix_templates{1,1:2};
  [POC_AL, POC_AR]= sediment_matrix_templates{1,1:2};
  [FeOOH_AL, FeOOH_AR]= sediment_matrix_templates{1,1:2};
  [FeS_AL, FeS_AR]= sediment_matrix_templates{1,1:2};
  [S8_AL, S8_AR]= sediment_matrix_templates{1,1:2};
  [FeS2_AL, FeS2_AR]= sediment_matrix_templates{1,1:2};
  [AlOH3_AL, AlOH3_AR]= sediment_matrix_templates{1,1:2};
  [Ca3PO42_AL, Ca3PO42_AR]= sediment_matrix_templates{1,1:2};
  [PO4adsa_AL, PO4adsa_AR]= sediment_matrix_templates{1,1:2};
  [PO4adsb_AL, PO4adsb_AR]= sediment_matrix_templates{1,1:2};
  [OMS_AL, OMS_AR]= sediment_matrix_templates{1,1:2};
  [FeOH3_AL, FeOH3_AR]= sediment_matrix_templates{1,1:2};
  [Chl_AL, Chl_AR]= sediment_matrix_templates{1,1:2};

  % Solute species:
  [O2_AL, O2_AR] = sediment_matrix_templates{2,1:2};
  [NO3_AL, NO3_AR] = sediment_matrix_templates{3,1:2};
  [SO4_AL, SO4_AR] = sediment_matrix_templates{4,1:2};
  [NH4_AL, NH4_AR] = sediment_matrix_templates{5,1:2};
  [Fe2_AL, Fe2_AR] = sediment_matrix_templates{6,1:2};
  [H2S_AL, H2S_AR] = sediment_matrix_templates{7,1:2};
  [S0_AL, S0_AR] = sediment_matrix_templates{8,1:2};
  [PO4_AL, PO4_AR] = sediment_matrix_templates{9,1:2};
  [Ca2_AL, Ca2_AR] = sediment_matrix_templates{10,1:2};
  [HS_AL, HS_AR] = sediment_matrix_templates{11,1:2};
  [H_AL, H_AR] = sediment_matrix_templates{12,1:2};
  [OH_AL, OH_AR] = sediment_matrix_templates{13,1:2};
  [CO2_AL, CO2_AR] = sediment_matrix_templates{14,1:2};
  [CO3_AL, CO3_AR] = sediment_matrix_templates{15,1:2};
  [HCO3_AL, HCO3_AR] = sediment_matrix_templates{16,1:2};
  [NH3_AL, NH3_AR] = sediment_matrix_templates{17,1:2};
  [H2CO3_AL, H2CO3_AR] = sediment_matrix_templates{18,1:2};
  [DOP_AL, DOP_AR] = sediment_matrix_templates{19,1:2};
  [DOC_AL, DOC_AR] = sediment_matrix_templates{20,1:2};
  [CH4_AL, CH4_AR] = sediment_matrix_templates{21,1:2};


  % Solving equations!!!
  % =========================================================================================================

  for i=2:m
    % =======================================================================================================
    % Solving Reaction eq-s
    % =======================================================================================================
    C0 = [O2(:,i-1), POP(:,i-1), POC(:,i-1), NO3(:,i-1), FeOH3(:,i-1), SO4(:,i-1), NH4(:,i-1), Fe2(:,i-1), FeOOH(:,i-1), H2S(:,i-1), HS(:,i-1), FeS(:,i-1), S0(:,i-1), PO4(:,i-1), S8(:,i-1), FeS2(:,i-1), AlOH3(:,i-1), PO4adsa(:,i-1), PO4adsb(:,i-1), Ca2(:,i-1), Ca3PO42(:,i-1), OMS(:,i-1), H(:,i-1), OH(:,i-1), CO2(:,i-1), CO3(:,i-1), HCO3(:,i-1), NH3(:,i-1), H2CO3(:,i-1), DOP(:,i-1), DOC(:,i-1), Chl(:,i-1), CH4(:,i-1)];

      if any(any(isnan(C0)))
          error('NaN')
      end



      int_method = 0;
      [C_new, rates(i-1)] = sediments_chemical_reactions_module(sediment_params, C0,dt, int_method);

      O2(:,i-1)      = C_new(:,1);
      POP(:,i-1)      = C_new(:,2);
      POC(:,i-1)     = C_new(:,3);
      NO3(:,i-1)     = C_new(:,4);
      FeOH3(:,i-1)   = C_new(:,5);
      SO4(:,i-1)     = C_new(:,6);
      NH4(:,i-1)     = C_new(:,7);
      Fe2(:,i-1)     = C_new(:,8);
      FeOOH(:,i-1)   = C_new(:,9);
      H2S(:,i-1)     = C_new(:,10);
      HS(:,i-1)      = C_new(:,11);
      FeS(:,i-1)     = C_new(:,12);
      S0(:,i-1)      = C_new(:,13);
      PO4(:,i-1)     = C_new(:,14);
      S8(:,i-1)      = C_new(:,15);
      FeS2(:,i-1)    = C_new(:,16);
      AlOH3(:,i-1)   = C_new(:,17);
      PO4adsa(:,i-1) = C_new(:,18);
      PO4adsb(:,i-1) = C_new(:,19);
      Ca2(:,i-1)     = C_new(:,20);
      Ca3PO42(:,i-1) = C_new(:,21);
      OMS(:,i-1)     = C_new(:,22);
      H(:,i-1)       = C_new(:,23);
      OH(:,i-1)      = C_new(:,24);
      CO2(:,i-1)     = C_new(:,25);
      CO3(:,i-1)     = C_new(:,26);
      HCO3(:,i-1)    = C_new(:,27);
      NH3(:,i-1)     = C_new(:,28);
      H2CO3(:,i-1)   = C_new(:,29);
      DOP(:,i-1)   = C_new(:,30);
      DOC(:,i-1)   = C_new(:,31);
      Chl(:,i-1)   = C_new(:,32);
      CH4(:,i-1)   = C_new(:,33);


    % =======================================================================================================
    % Solving Transport eq-s
    % =======================================================================================================

      O2(:,i) = pde_solver_dissolved(O2_AL, O2_AR, O2(:,i-1), sediment_bc.O2_c);
      POP(:,i) = pde_solver_solid(POP_AL, POP_AR, POP(:,i-1), sediment_bc.POP_fx, sediment_params.solid_flux_coef);
      POC(:,i) = pde_solver_solid(POC_AL, POC_AR, POC(:,i-1), sediment_bc.POC_fx, sediment_params.solid_flux_coef);
      NO3(:,i) = pde_solver_dissolved(NO3_AL, NO3_AR, NO3(:,i-1), sediment_bc.NO3_c);
      FeOH3(:,i) = pde_solver_solid(FeOH3_AL, FeOH3_AR, FeOH3(:,i-1), sediment_bc.FeOH3_fx, sediment_params.solid_flux_coef);
      SO4(:,i) = pde_solver_dissolved(SO4_AL, SO4_AR, SO4(:,i-1), sediment_bc.SO4_c);
      NH4(:,i) = pde_solver_dissolved(NH4_AL, NH4_AR, NH4(:,i-1), sediment_bc.NH4_c);
      Fe2(:,i) = pde_solver_dissolved(Fe2_AL, Fe2_AR, Fe2(:,i-1), sediment_bc.Fe2_c);
      FeOOH(:,i) = pde_solver_solid(FeOOH_AL, FeOOH_AR, FeOOH(:,i-1), sediment_bc.FeOOH_fx, sediment_params.solid_flux_coef);
      H2S(:,i) = pde_solver_dissolved(H2S_AL, H2S_AR, H2S(:,i-1), sediment_bc.H2S_c);
      HS(:,i) = pde_solver_dissolved(HS_AL, HS_AR, HS(:,i-1), sediment_bc.HS_c);
      FeS(:,i) = pde_solver_solid(FeS_AL, FeS_AR, FeS(:,i-1), sediment_bc.FeS_fx, sediment_params.solid_flux_coef);
      S0(:,i) = pde_solver_dissolved(S0_AL, S0_AR, S0(:,i-1), sediment_bc.S0_c);
      PO4(:,i) = pde_solver_dissolved(PO4_AL, PO4_AR, PO4(:,i-1), sediment_bc.PO4_c);
      S8(:,i) = pde_solver_solid(S8_AL, S8_AR, S8(:,i-1), sediment_bc.S8_fx, sediment_params.solid_flux_coef);
      FeS2(:,i) = pde_solver_solid(FeS2_AL, FeS2_AR, FeS2(:,i-1), sediment_bc.FeS2_fx, sediment_params.solid_flux_coef);
      AlOH3(:,i) = pde_solver_solid(AlOH3_AL, AlOH3_AR, AlOH3(:,i-1), sediment_bc.AlOH3_fx, sediment_params.solid_flux_coef);
      PO4adsa(:,i) = pde_solver_solid(PO4adsa_AL, PO4adsa_AR, PO4adsa(:,i-1), sediment_bc.PO4adsa_fx, sediment_params.solid_flux_coef);
      PO4adsb(:,i) = pde_solver_solid(PO4adsb_AL, PO4adsb_AR, PO4adsb(:,i-1), sediment_bc.PO4adsb_fx, sediment_params.solid_flux_coef);
      Ca2(:,i) = pde_solver_dissolved(Ca2_AL, Ca2_AR, Ca2(:,i-1), sediment_bc.Ca2_c);
      Ca3PO42(:,i) = pde_solver_solid(Ca3PO42_AL, Ca3PO42_AR, Ca3PO42(:,i-1), sediment_bc.Ca3PO42_fx, sediment_params.solid_flux_coef);
      OMS(:,i) = pde_solver_solid(OMS_AL, OMS_AR, OMS(:,i-1), sediment_bc.OMS_fx, sediment_params.solid_flux_coef);
      H(:,i) = pde_solver_dissolved(H_AL, H_AR, H(:,i-1), sediment_bc.H_c);
      OH(:,i) = pde_solver_dissolved(OH_AL, OH_AR, OH(:,i-1), sediment_bc.OH_c);
      CO2(:,i) = pde_solver_dissolved(CO2_AL, CO2_AR, CO2(:,i-1), sediment_bc.CO2_c);
      CO3(:,i) = pde_solver_dissolved(CO3_AL, CO3_AR, CO3(:,i-1), sediment_bc.CO3_c);
      HCO3(:,i) = pde_solver_dissolved(HCO3_AL, HCO3_AR, HCO3(:,i-1), sediment_bc.HCO3_c);
      NH3(:,i) = pde_solver_dissolved(NH3_AL, NH3_AR, NH3(:,i-1), sediment_bc.NH3_c);
      H2CO3(:,i) = pde_solver_dissolved(H2CO3_AL, H2CO3_AR, H2CO3(:,i-1), sediment_bc.H2CO3_c);
      DOP(:,i) = pde_solver_dissolved(DOP_AL, DOP_AR, DOP(:,i-1), sediment_bc.DOP_c);
      DOC(:,i) = pde_solver_dissolved(DOC_AL, DOC_AR, DOC(:,i-1), sediment_bc.DOC_c);
      Chl(:,i) = pde_solver_solid(Chl_AL, Chl_AR, Chl(:,i-1), sediment_bc.Chl_fx, sediment_params.solid_flux_coef);
      CH4(:,i) = pde_solver_dissolved(CH4_AL, CH4_AR, CH4(:,i-1), sediment_bc.CH4_c);

      % Estimate fluxes:

      sediment_bioirrigation_fluxes.O2(i-1)   = integrate_over_depth_2( bioirrigation(O2(:, i), alfax, phi), x);
      sediment_bioirrigation_fluxes.PO4(i-1)  = integrate_over_depth_2( bioirrigation(PO4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.Fe2(i-1)  = integrate_over_depth_2( bioirrigation(Fe2(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.NO3(i-1)  = integrate_over_depth_2( bioirrigation(NO3(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.NH4(i-1)  = integrate_over_depth_2( bioirrigation(NH4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.SO4(i-1)  = integrate_over_depth_2( bioirrigation(SO4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.DOP(i-1) = integrate_over_depth_2( bioirrigation(DOP(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.DOC(i-1) = integrate_over_depth_2( bioirrigation(DOC(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CH4(i-1) = integrate_over_depth_2( bioirrigation(CH4(:, i),  alfax,  phi), x);

      sediment_transport_fluxes.POP(i-1)          = -sediment_bc.POP_fx;
      sediment_transport_fluxes.Chl(i-1)          = -sediment_bc.Chl_fx;
      sediment_transport_fluxes.POC(i-1)          = -sediment_bc.POC_fx;
      sediment_transport_fluxes.FeOH3(i-1)        = -sediment_bc.FeOH3_fx;
      sediment_transport_fluxes.AlOH3(i-1)        = -sediment_bc.AlOH3_fx;
      sediment_transport_fluxes.PO4adsa(i-1)      = -sediment_bc.PO4adsa_fx;
      sediment_transport_fluxes.PO4adsb(i-1)      = -sediment_bc.PO4adsb_fx;
      sediment_transport_fluxes.O2(i-1)           = top_sediment_diffusion_flux(O2(:, i), D_O2, dx, phi);
      sediment_transport_fluxes.PO4(i-1)          = top_sediment_diffusion_flux(PO4(:, i), D_PO4, dx, phi);
      sediment_transport_fluxes.NO3(i-1)          = top_sediment_diffusion_flux(NO3(:, i), D_NO3, dx, phi);
      sediment_transport_fluxes.Fe2(i-1)          = top_sediment_diffusion_flux(Fe2(:, i), D_Fe2, dx, phi);
      sediment_transport_fluxes.NH4(i-1)          = top_sediment_diffusion_flux(NH4(:, i), D_NH4, dx, phi);
      sediment_transport_fluxes.SO4(i-1)          = top_sediment_diffusion_flux(SO4(:, i), D_SO4, dx, phi);
      sediment_transport_fluxes.DOP(i-1)         = top_sediment_diffusion_flux(DOP(:, i), D_DOP, dx, phi);
      sediment_transport_fluxes.DOC(i-1)         = top_sediment_diffusion_flux(DOC(:, i), D_DOC, dx, phi);
      sediment_transport_fluxes.CH4(i-1)         = top_sediment_diffusion_flux(CH4(:, i), D_DOC, dx, phi);


    % pH Module
    if sediment_params.pH_algorithm ~= 0
      [H(:,i), OH(:,i), DOC(:,i), HCO3(:,i), CO2(:,i), CO3(:,i), NH3(:,i), NH4(:,i), HS(:,i), H2S(:,i)] = pH_module(sediment_params.pH_algorithm, H(:,i), OH(:,i), H2CO3(:,i), HCO3(:,i), CO2(:,i), CO3(:,i), NH3(:,i), NH4(:,i), HS(:,i), H2S(:,i), Fe2(:,i), Ca2(:,i), NO3(:,i), SO4(:,i), PO4(:,i), FeS(:,i), FeS2(:,i), FeOH3(:,i), FeOOH(:,i), Ca3PO42(:,i), PO4adsa(:,i), PO4adsb(:,i), sediment_bc.T);
    end

  end

% Estimate flux
  sediment_transport_fluxes.O2           = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.O2), 31998);
  sediment_transport_fluxes.POP          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.POP), 30973.762);
  sediment_transport_fluxes.POC          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.POC), 12010.7);
  sediment_transport_fluxes.PO4          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.PO4), 30973.762);
  sediment_transport_fluxes.NO3          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.NO3), 62004);
  sediment_transport_fluxes.FeOH3        = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.FeOH3), 106867.0);
  sediment_transport_fluxes.Fe2          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.Fe2), 55845);
  sediment_transport_fluxes.NH4          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.NH4), 18038);
  sediment_transport_fluxes.AlOH3        = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.AlOH3), 78003.6);
  sediment_transport_fluxes.PO4adsa      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.PO4adsa), 30973.762);
  sediment_transport_fluxes.PO4adsb      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.PO4adsb), 30973.762);
  sediment_transport_fluxes.SO4          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.SO4), 96062);
  sediment_transport_fluxes.DOP         = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.DOP), 30973.762);
  sediment_transport_fluxes.DOC         = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.DOC), 12010.7);
  sediment_transport_fluxes.Chl          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.Chl), 30973.762);

  sediment_bioirrigation_fluxes.O2   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.O2), 31998);
  sediment_bioirrigation_fluxes.PO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.PO4), 30973.762);
  sediment_bioirrigation_fluxes.Fe2  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.Fe2), 55845);
  sediment_bioirrigation_fluxes.NO3  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.NO3), 62004);
  sediment_bioirrigation_fluxes.NH4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.NH4), 18038);
  sediment_bioirrigation_fluxes.SO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.SO4), 96062);
  sediment_bioirrigation_fluxes.DOP = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.DOP), 30973.762);
  sediment_bioirrigation_fluxes.DOC = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.DOC), 12010.7);
  sediment_bioirrigation_fluxes.CH4 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CH4), 16042.5);


  sediment_concentrations.O2 = O2(:,end);
  sediment_concentrations.POP = POP(:,end);
  sediment_concentrations.POC = POC(:,end);
  sediment_concentrations.NO3 = NO3(:,end);
  sediment_concentrations.FeOH3 = FeOH3(:,end);
  sediment_concentrations.SO4 = SO4(:,end);
  sediment_concentrations.NH4 = NH4(:,end);
  sediment_concentrations.Fe2 = Fe2(:,end);
  sediment_concentrations.FeOOH = FeOOH(:,end);
  sediment_concentrations.H2S = H2S(:,end);
  sediment_concentrations.HS = HS(:,end);
  sediment_concentrations.FeS = FeS(:,end);
  sediment_concentrations.S0 = S0(:,end);
  sediment_concentrations.PO4 = PO4(:,end);
  sediment_concentrations.S8 = S8(:,end);
  sediment_concentrations.FeS2 = FeS2(:,end);
  sediment_concentrations.AlOH3 = AlOH3(:,end);
  sediment_concentrations.PO4adsa = PO4adsa(:,end);
  sediment_concentrations.PO4adsb = PO4adsb(:,end);
  sediment_concentrations.Ca2 = Ca2(:,end);
  sediment_concentrations.Ca3PO42 = Ca3PO42(:,end);
  sediment_concentrations.OMS = OMS(:,end);
  sediment_concentrations.H = H(:,end);
  sediment_concentrations.OH = OH(:,end);
  sediment_concentrations.CO2 = CO2(:,end);
  sediment_concentrations.CO3 = CO3(:,end);
  sediment_concentrations.HCO3 = HCO3(:,end);
  sediment_concentrations.NH3 = NH3(:,end);
  sediment_concentrations.H2CO3 = H2CO3(:,end);
  sediment_concentrations.DOP = DOP(:,end);
  sediment_concentrations.DOC = DOC(:,end);
  sediment_concentrations.Chl = Chl(:,end);
  sediment_concentrations.CH4 = CH4(:,end);

  sediment_concentrations.pH = -log10(H(:,end)*10^-3);


  % Estimate average rate during the day
  if sediment_params.rate_estimator_switch
    fields = fieldnames(rates);
    for i = 1:numel(fields)
        r.(fields{i}) = 0;
        for j=1:m-1
            r.(fields{i}) = r.(fields{i}) + rates(j).(fields{i});
        end
        r.(fields{i}) = r.(fields{i})/(m-1);
    end

    sediment_additional_results.rates = r;

  else
    sediment_additional_results.rates = false;
  end

    if any(isnan(sediment_transport_fluxes.O2))| any(isnan(sediment_bc.POP_fx))| any(isnan(sediment_bc.POC_fx))| any(isnan(sediment_bc.FeOH3_fx))| any(isnan(O2)) | any(isnan(POP)) | any(isnan(POC)) | any(isnan(NO3)) | any(isnan(FeOH3)) | any(isnan(SO4)) | any(isnan(NH4)) | any(isnan(Fe2)) | any(isnan(FeOOH)) | any(isnan(H2S)) | any(isnan(HS)) | any(isnan(FeS)) | any(isnan(S0)) | any(isnan(PO4)) | any(isnan(S8)) | any(isnan(FeS2)) | any(isnan(AlOH3)) | any(isnan(PO4adsa)) | any(isnan(PO4adsb)) | any(isnan(H)) | any(isnan(Ca2)) | any(isnan(Ca3PO42)) | any(isnan(OMS)) | any(isnan(OH)) | any(isnan(HCO3)) | any(isnan(CO2)) | any(isnan(CO3)) | any(isnan(NH3)) | any(isnan(H2CO3)) | any(isnan(Chl))  any(isnan(CH4))
      error('Breaking out of Sediments function: NaN values');
    end


end



function [H, OH, H2CO3, HCO3, CO2, CO3, NH3, NH4, HS, H2S] = pH_module(algorithm, H, OH, H2CO3, HCO3, CO2, CO3, NH3, NH4, HS, H2S, Fe2, Ca2, NO3, SO4, PO4, FeS, FeS2, FeOH3, FeOOH, Ca3PO42, PO4adsa, PO4adsb, Temperature)
  %% pH_module: pH equilibrium function
  % 0. No pH module
  % 1. Stumm, W. & Morgan, J., 1995. Aquatic Chemistry. implemented in MATLAB
  % 2. Stumm, W. & Morgan, J., 1995. Aquatic Chemistry. implemented in C++
  % 3. Phreeqc
  % 4. Delta function (under construction)

  % NOTE: First point is boundary condition therefore start FOR loop from 2:end

    if algorithm == 1 % Stumm, W. & Morgan, J., 1995. Aquatic Chemistry. implemented in MATLAB
      % The fasted algorithm is Levenberg-Marquardt
      options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','TolFun',1e-7, 'TolX',1e-7); % Option to display output

      for i=1:size(H,1)   % Stumm, W. & Morgan, J., 1995. Aquatic Chemistry. implemented in MATLAB
        % initial guess
        x = [sqrt(H(i)); sqrt(HCO3(i)); sqrt(CO2(i)); sqrt(CO3(i)); sqrt(NH3(i)); sqrt(NH4(i)); sqrt(HS(i)); sqrt(H2S(i)); sqrt(OH(i)); sqrt(H2CO3(i))];
        [x,~] = fsolve(@SM_eqs, x, options, sqrt(H(i)), sqrt(OH(i)), sqrt(H2CO3(i)), sqrt(HCO3(i)), sqrt(CO2(i)), sqrt(CO3(i)), sqrt(NH3(i)), sqrt(NH4(i)), sqrt(HS(i)), sqrt(H2S(i)), sqrt(Fe2(i)), sqrt(Ca2(i)), sqrt(NO3(i)), sqrt(SO4(i)), sqrt(PO4(i))); % Call solver
        H(i) = x(1)^2;
        HCO3(i) = x(2)^2;
        CO2(i) = x(3)^2;
        CO3(i) = x(4)^2;
        NH3(i) = x(5)^2;
        NH4(i) = x(6)^2;
        HS(i) = x(7)^2;
        H2S(i) = x(8)^2;
        OH(i) = x(9)^2;
        H2CO3(i) = x(10)^2;
      end


    elseif algorithm == 2 % Stumm, W. & Morgan, J., 1995. Aquatic Chemistry. implemented in C++
      for i=1:size(H,1)
        % NOTE: We have NaN if start from 1 (with BC).
        in =[H(i) HCO3(i) CO2(i) CO3(i) NH3(i) NH4(i) HS(i) H2S(i) OH(i) H2CO3(i) Fe2(i) Ca2(i) NO3(i) SO4(i) PO4(i)];
        [out] = pH(in);
        H(i) = out(1);
        HCO3(i) = out(2);
        CO2(i) = out(3);
        CO3(i) = out(4);
        NH3(i) = out(5);
        NH4(i) = out(6);
        HS(i) = out(7);
        H2S(i) = out(8);
        OH(i) = out(9);
        H2CO3(i) = out(10);
      end

    elseif algorithm == 3 %
        in =[H HCO3 CO2 CO3 NH3 NH4 HS H2S OH H2CO3 Fe2 Ca2 NO3 SO4 PO4 FeS FeS2 FeOH3 FeOOH Ca3PO42 PO4adsa PO4adsb];
        [pH_est] = pH_phreeqc(size(H,1),in);

        Kc1=5.01*10^(-7); Kc2=4.78*10^(-11); Knh=5.62*10^(-10); Khs=1.3*10^(-7); Kw=10^(-14); Kc0 = 1.7*10^(-3);
        H = 10.^(-pH_est');
        Ct = H2CO3 + HCO3 + CO3;
        Nt = NH3 + NH4;
        St = HS + H2S;
        OH = Kw./H;

        % TODO: Carbonate equilibrium doesn't work now. You can  estimate fractions based on pH after model run
        % Tz = Temperature + 273.15;
        % K0 = -60.2409+93.4517.*(100./Tz)+23.3585*log(Tz/100); %~mol/(kg*atm)
        % K0 = exp(K0); %CO2; mol/(kg*atm)
        % K1 = 290.9097-14554.21./Tz-45.0575*log(Tz); %~mol/kg
        % K1 = exp(K1); %HCO3; mol/kg
        % K2 = 207.6548-11843.79./Tz-33.6485*log(Tz); %~mol/kg
        % K2 = exp(K2); %CO3; mol/kg
        % Kw = 148.9802-13847.26./Tz-23.6521*log(Tz); %~mol/kg
        % Kw = exp(Kw); %H2O; (mol/kg)^2
        % CO2mfrac = H.*H./((H.*H+H.*K1+K1.*K2)); %mol CO2 / mol DIC
        % HCO3mfrac = H.*K1./((H.*H+H.*K1+K1.*K2)); %-mol HCO3 / mol DIC
        % CO3mfrac = K1.*K2./((H.*H+H.*K1+K1.*K2)); %-mol CO3 / mol DIC
        % H2CO3 = CO2mfrac .* Ct;
        % HCO3 = HCO3mfrac .* Ct;
        % CO3 = CO3mfrac .* Ct;
        % H2CO3 = (1 + Kc1./H + Kc1*Kc2./H.^2).^-1 .* Ct;
        % CO2 = H2CO3./Kc0; % CO2(aq)
        % HCO3 = (H./Kc1 + 1 + Kc2./H).^-1 .* Ct;
        % CO3 = (H.^2./Kc1./Kc2 + H./Kc2 + 1).^-1 .* Ct;

        NH4 = (H./ (Knh + H)) .* Nt;
        NH3 = Nt - NH4;
        H2S = (H./ (Khs + H)) .* St;
        HS = St - H2S;
        H = H*1e3;
        OH = OH*1e3;



    elseif algorithm == 4 % Delta function
      options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','TolFun',1e-7, 'TolX',1e-7); % Option to display output
      delta = [0,0,0,0,0,0,0,0,0];
      for i=1:size(H,1)
        [delta,fval] = fsolve(@delta_eqs, delta, options, H(i), OH(i), HCO3(i), CO2(i), CO3(i), NH3(i), NH4(i), HS(i), H2S(i)); % Call solver
        H(i) = H(i) + delta(1)+delta(2)+delta(3)+delta(4)+delta(5);
        OH(i) = OH(i) + delta(5);
        HCO3(i) = HCO3(i) + delta(1) - delta(2);
        CO3(i) = CO3(i) + delta(2);
        CO2(i) = CO2(i) - delta(1);
        NH3(i) = NH3(i) + delta(3);
        NH4(i) = NH4(i) - delta(3);
        HS(i) = HS(i) + delta(4);
        H2S(i) = H2S(i) - delta(4);
      end
    end

end


function F = SM_eqs(x, H, OH, H2CO3, HCO3, CO2, CO3, NH3, NH4, HS, H2S, Fe2, Ca2, NO3, SO4, PO4)
  %% Equations according to Stumm, W. & Morgan, J., 1995. Aquatic Chemistry implemented in MATLAB.

  % x(1)^2=H;   % x(2)^2=HCO3;  % x(3)^2=CO2;  % x(4)^2=CO3;  % x(5)^2=NH3;  % x(6)^2=NH4;  % x(7)^2=HS;  % x(8)^2=H2S;  % x(9)^2=OH;

  Kc1=5.01*10^(-7).*10^3; Kc2=4.78*10^(-11).*10^3; Knh=5.62*10^(-10).*10^3; Khs=1.3*10^(-7).*10^3; Kw=10^(-14).*10^6; Kc0 = 1.7*10^(-3);

  F = [
      x(10)^2 - Kc0*x(3)^2;
      x(1)^2 * x(2)^2 - Kc1*x(10)^2;
      x(1)^2 * x(4)^2 - Kc2*x(2)^2;
      x(1)^2 * x(5)^2 - Knh*x(6)^2;
      x(1)^2 * x(7)^2 - Khs*x(8)^2;
      x(1)^2 * x(9)^2 - Kw;
     % mass balance
     x(5)^2 + x(6)^2 - NH4^2 - NH3^2;
     x(7)^2 + x(8)^2 - HS^2 - H2S^2;
     x(4)^2 + x(2)^2 + x(3)^2 + x(10)^2 - CO3^2 - HCO3^2 - CO2^2 - H2CO3^2;
     % charge balance
     x(1)^2 + x(6)^2 + 2*Fe2^2 + 2*Ca2^2 - (x(2)^2 + 2*x(4)^2 + x(7)^2 + x(9)^2 + NO3^2 + 2*SO4^2 + 3*PO4^2);
     ];
end

function F = delta_eqs(delta, H, OH, HCO3, CO2, CO3, NH3, NH4, HS, H2S)
  %% equations according to Delta function
  Kc1=10^(-6.4).*10^3; Kc2=10^(-10.3).*10^3; Knh=10^(-9.3).*10^3; Khs=10^(-7).*10^6; Kw=10^(-14).*10^6;

  F = [(H + delta(1) + delta(2) + delta(3) + delta(4) + delta(5) )* (HCO3 + delta(1) - delta(2) ) - Kc1*(CO2 - delta(1));
     (H + delta(1) + delta(2) + delta(3) + delta(4) + delta(5) )* (CO3 + delta(2) ) -  Kc2*(HCO3 - delta(2) + delta(1));
     (H + delta(1) + delta(2) + delta(3) + delta(4) + delta(5) )* (NH3 + delta(3) ) - Knh* (NH4 - delta(3));
     (H + delta(1) + delta(2) + delta(3) + delta(4) + delta(5) )* (HS + delta(4) ) - Khs* (H2S - delta(4));
     (H + delta(1) + delta(2) + delta(3) + delta(4) + delta(5) )* (OH + delta(5) ) - Kw;
     ];
end

function C_new = pde_solver_dissolved(AL, AR, C_old, const_bc)
    C_old(1) = const_bc;
    temp = AR*C_old;
    % temp(1) = const_bc;
    C_new = AL\ temp;
    % C_new(1) = const_bc;
    C_new = (C_new>0).*C_new;
end


function C_new = pde_solver_solid(AL, AR, C_old, flux_bc, coef)
      temp = AR*C_old;
      temp(1) = temp(1) + flux_bc * coef;
      C_new = AL\ temp;
      C_new = (C_new>0).*C_new;
end

function [C_new, rates] = sediments_chemical_reactions_module(sediment_params, C0,dt, method)
    if method == 0
        [C_new, rates] = rk4(sediment_params, C0,dt);
    elseif method == 1
        [C_new, rates] = butcher5(sediment_params, C0,dt);
    end
    C_new = (C_new>0).*C_new;
end

%% rk4: Runge-Kutta 4th order integration
function [C_new, rates] = rk4(sediment_params, C0, dt)
    % ts - how many time steps during 1 day

    [dcdt_1, r_1] = sediment_rates(sediment_params, C0, dt);
    k_1 = dt.*dcdt_1;
    [dcdt_2, r_2] = sediment_rates(sediment_params, C0+0.5.*k_1, dt);
    k_2 = dt.*dcdt_2;
    [dcdt_3, r_3] = sediment_rates(sediment_params, C0+0.5.*k_2, dt);
    k_3 = dt.*dcdt_3;
    [dcdt_4, r_4] = sediment_rates(sediment_params, C0+k_3, dt);
    k_4 = dt.*dcdt_4;
    C_new = C0 + (k_1+2.*k_2+2.*k_3+k_4)/6;
    C0 = C_new;

    if sediment_params.rate_estimator_switch
      % average rate
      fields = fieldnames(r_1);
      for fld_idx = 1:numel(fields)
        rates.(fields{fld_idx}) = (r_1.(fields{fld_idx}) + 2*r_2.(fields{fld_idx}) + 2*r_3.(fields{fld_idx}) + r_4.(fields{fld_idx}))/6;
      end
    else
      rates = false;
    end
end

%% butcher5: Butcher's Fifth-Order Runge-Kutta
function [C_new, rates] = butcher5(sediment_params, C0,dt)

    [dcdt_1, r_1] = sediment_rates(sediment_params, C0, dt);
    k_1 = dt.*dcdt_1;
    [dcdt_2, r_2] = sediment_rates(sediment_params, C0 + 1/4.*k_1, dt);
    k_2 = dt.*dcdt_2;
    [dcdt_3, r_3] = sediment_rates(sediment_params, C0 + 1/8.*k_1 + 1/8.*k_2, dt);
    k_3 = dt.*dcdt_3;
    [dcdt_4, r_4] = sediment_rates(sediment_params, C0 - 1/2.*k_2 + k_3, dt);
    k_4 = dt.*dcdt_4;
    [dcdt_5, r_5] = sediment_rates(sediment_params, C0 + 3/16.*k_1 + 9/16.*k_4, dt);
    k_5 = dt.*dcdt_5;
    [dcdt_6, r_6] = sediment_rates(sediment_params, C0 - 3/7.*k_1 + 2/7.*k_2 + 12/7.*k_3 - 12/7.*k_4 + 8/7.*k_5, dt);
    k_6 = dt.*dcdt_6;
    C_new = C0 + (7.*k_1 + 32.*k_3 + 12.*k_4 + 32.*k_5 + 7.*k_6)/90;
    C0 = C_new;

    % average rate
    if sediment_params.rate_estimator_switch
      fields = fieldnames(r_1);
      for fld_idx = 1:numel(fields)
        rates.(fields{fld_idx}) = (7*r_1.(fields{fld_idx}) + 32*r_3.(fields{fld_idx}) + 12*r_4.(fields{fld_idx}) + 32*r_5.(fields{fld_idx}) + 7*r_6.(fields{fld_idx}))/90;
      end
    else
      rates = false;
    end
end

%% top_sediment_rate_to_flux: returns the flux of species at SWI converted to th units used in WC [ mg m-2 d-1 ].
function [flux] = top_sediment_rate_to_flux(R, dx)
    % TODO: Check units here!!
    % R - rate [umol cm-3 y-1]
    % dx - mesh size [cm];
    flux = integrate_over_depth(R, dx);
end

function [int_rate] = integrate_over_depth(R, dx)
  %% integrate_over_depth: integrates the rates of reaction over the depth and return the average values for the current day
  % int_rate  [umol cm-2 year -1 ]
  % R - rate of interest
  % dx - the mesh size
  % m - number of time step in 1 day
  int_rate = sum(daily_average(R),1)*dx;
end

function [int_rate] = integrate_over_depth_2(R, z)
  %% integrate_over_depth_2: integrates the rates of reaction over the depth and return the average values for the current day using trapezoidal rule
  % int_rate  [umol cm-2 year -1 ]
  % R - rate of interest
  % z - the depth
  if size(R,1) == 1
    int_rate = 0;
  else
    int_rate = trapz(z,R);
  end
end

%% daily_average: returns the average rate during 1 run (usually it is during 1 day if run with MyLake)
function [averaged] = daily_average(R)
  % R - rate of interest
  averaged = sum(R,2)/size(R,2);
end

function [flux] = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(fx, M_C)
  flux = fx * M_C * 10^4 / 365 / 10^6; % [umol/cm^2/y] -> [mg/m2/d]
end

function [flux] = top_sediment_diffusion_flux(C, D, dx, phi)
  % calculates flux of the particular dissolved specie through the top boundary of the sediment
  % in [ mg m-2 d-1 ] units
  % C(1) - BC of dissolved species
  % C - concentration
  % D - diffusion coefficient
  % M_C - molar mass in [ mg mol-1]
  % phi - porosity (no porosity because C is the concentration in pores (not bulk))

  % fourth-order
  flux = D * (-25 * phi(2)*C(2) + 48 * phi(3)*C(3) - 36 * phi(4)*C(4) + 16 * phi(5)*C(5) - 3 * phi(6)*C(6)) / dx / 12;  %  [umol/cm^2/y]

  % third order
  % flux = D * (-11 * C(1) + 18 * C(2) - 9 * C(3) + 2 * C(4)) / dx / 6;  %  [umol/cm^2/y]
  % flux = 0;  %  [umol/cm^2/y]

  % second order
  % flux = D * (-3 * C(1) + 4 * C(2) - C(3)) / dx / 2;  %  [umol/cm^2/y]

  % first order
  % flux = D * (C(3) - C(1)) / 2 / dx;  %  [umol/cm^2/y]
end

function bioR = bioirrigation(C, alfax, phi)
  % bioirrigation rate is the "artificial" function represents the bioirrigation by organisms in the sediment (worms etc) implemented according to Boudreau, B.P., 1999.
  % Co - bc wc-sediment value of current species
  % C - concentration profile of current species
  % phi - porosity
  Co = C(1);
  bioR = alfax .* (C - Co); % .* fi
  % NOTE:Disabled?
  % bioR = 0;
end


function [dcdt, r] = sediment_rates(sediment_params, C, dt)
% parameters for water-column chemistry
% NOTE: the rates are the same as in sediments. Units are per "year" due to time step is in year units too;


    dcdt=zeros(size(C));

    if any(isnan(C))
      error('Breaking out of Sediments function: NaN values');
    end

    Ox      = C(:,1) .* (C(:,1)>0) ;
    POP      = C(:,2) .* (C(:,2)>0) ;
    POC     = C(:,3) .* (C(:,3)>0) ;
    NO3     = C(:,4) .* (C(:,4)>0) ;
    FeOH3   = C(:,5) .* (C(:,5)>0) ;
    SO4     = C(:,6) .* (C(:,6)>0) ;
    NH4     = C(:,7) .* (C(:,7)>0) ;
    Fe2     = C(:,8) .* (C(:,8)>0) ;
    FeOOH   = C(:,9) .* (C(:,9)>0) ;
    H2S     = C(:,10) .* (C(:,10)>0) ;
    HS      = C(:,11) .* (C(:,11)>0) ;
    FeS     = C(:,12) .* (C(:,12)>0) ;
    S0      = C(:,13) .* (C(:,13)>0) ;
    PO4     = C(:,14) .* (C(:,14)>0) ;
    S8      = C(:,15) .* (C(:,15)>0) ;
    FeS2    = C(:,16) .* (C(:,16)>0) ;
    AlOH3   = C(:,17) .* (C(:,17)>0) ;
    PO4adsa = C(:,18) .* (C(:,18)>0) ;
    PO4adsb = C(:,19) .* (C(:,19)>0) ;
    Ca2     = C(:,20) .* (C(:,20)>0) ;
    Ca3PO42 = C(:,21) .* (C(:,21)>0) ;
    OMS     = C(:,22) .* (C(:,22)>0) ;
    H       = C(:,23) .* (C(:,23)>0) ;
    OH      = C(:,24) .* (C(:,24)>0) ;
    CO2     = C(:,25) .* (C(:,25)>0) ;
    CO3     = C(:,26) .* (C(:,26)>0) ;
    HCO3    = C(:,27) .* (C(:,27)>0) ;
    NH3     = C(:,28) .* (C(:,28)>0) ;
    H2CO3   = C(:,29) .* (C(:,29)>0) ;
    DOP    = C(:,30) .* (C(:,30)>0) ;
    DOC    = C(:,31) .* (C(:,31)>0) ;
    Chl    = C(:,32) .* (C(:,32)>0) ;
    CH4    = C(:,33) .* (C(:,33)>0) ;

    k_Chl =  sediment_params.k_Chl;
    k_POP =  sediment_params.k_POP;
    k_POC = sediment_params.k_POC;
    k_DOP = sediment_params.k_DOP;
    k_DOC = sediment_params.k_DOC;
    Km_O2 = sediment_params.Km_O2;
    Km_NO3 = sediment_params.Km_NO3;
    Km_FeOH3 = sediment_params.Km_FeOH3;
    Km_FeOOH = sediment_params.Km_FeOOH;
    Km_SO4 = sediment_params.Km_SO4;
    Km_oxao = sediment_params.Km_oxao;
    Km_amao = sediment_params.Km_amao;
    Kin_O2 = sediment_params.Kin_O2;
    Kin_NO3  = sediment_params.Kin_NO3;
    Kin_FeOH3 = sediment_params.Kin_FeOH3;
    Kin_FeOOH = sediment_params.Kin_FeOOH;
    k_amox = sediment_params.k_amox;
    k_Feox = sediment_params.k_Feox;
    k_Sdis = sediment_params.k_Sdis;
    k_Spre = sediment_params.k_Spre;
    k_FeS2pre = sediment_params.k_FeS2pre;
    k_pdesorb_c = sediment_params.k_pdesorb_c;
    k_pdesorb_a = sediment_params.k_pdesorb_a;
    k_pdesorb_b = sediment_params.k_pdesorb_b;
    % k_alum = sediment_params.k_alum;
    k_rhom   = sediment_params.k_rhom;
    k_tS_Fe = sediment_params.k_tS_Fe;
    Ks_FeS = sediment_params.Ks_FeS;
    k_Fe_dis = sediment_params.k_Fe_dis;
    k_Fe_pre = sediment_params.k_Fe_pre;
    k_apa  = sediment_params.k_apa;
    kapa = sediment_params.kapa;
    k_oms = sediment_params.k_oms;
    k_tsox = sediment_params.k_tsox;
    k_FeSpre = sediment_params.k_FeSpre;
    f_pfe = sediment_params.f_pfe;
    accel = sediment_params.accel;
    Cx1   = sediment_params.Cx1;
    Ny1   = sediment_params.Ny1;
    Pz1   = sediment_params.Pz1;
    Cx2   = sediment_params.Cx2;
    Ny2   = sediment_params.Ny2;
    Pz2   = sediment_params.Pz2;
    Cx3   = sediment_params.Cx3;
    Ny3   = sediment_params.Ny3;
    Pz3   = sediment_params.Pz3;
    alfax = sediment_params.alfax;
    phi    = sediment_params.phi;

    tot_FeOH3 = PO4adsa + FeOH3;

    f_O2    = Ox ./  (Km_O2 + Ox);
    f_NO3   = NO3 ./  (Km_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    f_FeOH3 = tot_FeOH3 ./  (Km_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    f_FeOOH = FeOOH ./  (Km_FeOOH + FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    f_SO4   = SO4 ./ (Km_SO4 + SO4 ) .* Kin_FeOOH ./ (Kin_FeOOH + FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    f_CH4 = 1 - f_O2 - f_NO3 - f_FeOH3 - f_FeOOH - f_SO4;

    Sum_H2S = H2S + HS;;

    Sat_FeS = Fe2*1e-3 .* Sum_H2S*1e-3 ./ (H*1e-3+1e-16).^2 ./ Ks_FeS;

    part_PO4ads_tot_Fe_a = PO4adsa ./ (tot_FeOH3+1e-16); % ratio of ads P to total Fe(III)

    R1a = k_POP.*POP .* f_O2 * accel;
    R1b = k_POC.*POC .* f_O2 * accel;
    R1c = k_DOP.* DOP .* f_O2 * accel;
    R1d = k_DOC.*DOC .* f_O2 * accel;
    R1f = k_Chl.*Chl .* f_O2 * accel;

    R2a = k_POP.*POP .* f_NO3 * accel;
    R2b = k_POC.*POC .* f_NO3 * accel;
    R2c = k_DOP.*DOP .* f_NO3 * accel;
    R2d = k_DOC.*DOC .* f_NO3 * accel;
    R2f = k_Chl.*Chl .* f_NO3 * accel;

    R3a_Fe = (1 - part_PO4ads_tot_Fe_a) .* k_POP.*POP .* f_FeOH3;
    R3b_Fe = (1 - part_PO4ads_tot_Fe_a) .* k_POC .*POC .* f_FeOH3;
    R3c_Fe = (1 - part_PO4ads_tot_Fe_a) .* k_DOP .*DOP .* f_FeOH3;
    R3d_Fe = (1 - part_PO4ads_tot_Fe_a) .* k_DOC .*DOC .* f_FeOH3;
    R3f_Fe = (1 - part_PO4ads_tot_Fe_a) .* k_Chl .*Chl .* f_FeOH3;
    R3a_P = part_PO4ads_tot_Fe_a .* k_POP.*POP .* f_FeOH3;
    R3b_P = part_PO4ads_tot_Fe_a .* k_POC .*POC .* f_FeOH3;
    R3c_P = part_PO4ads_tot_Fe_a .* k_DOP .*DOP .* f_FeOH3;
    R3d_P = part_PO4ads_tot_Fe_a .* k_DOC .*DOC .* f_FeOH3;
    R3f_P = part_PO4ads_tot_Fe_a .* k_Chl .*Chl .* f_FeOH3;
    R3a = R3a_Fe + R3a_P;
    R3b = R3b_Fe + R3b_P;
    R3c = R3c_Fe + R3c_P;
    R3d = R3d_Fe + R3d_P;
    R3f = R3f_Fe + R3f_P;

    R4a = k_POP.*POP .* f_FeOOH;
    R4b = k_POC.*POC .* f_FeOOH;
    R4c = k_DOP.*DOP .* f_FeOOH;
    R4d = k_DOC.*DOC .* f_FeOOH;
    R4f = k_Chl.*Chl .* f_FeOOH;

    R5a = k_POP.*POP .* f_SO4  ;
    R5b = k_POC.*POC .* f_SO4 ;
    R5c = k_DOP.*DOP .* f_SO4 ;
    R5d = k_DOC.*DOC .* f_SO4 ;
    R5f = k_Chl.*Chl .* f_SO4 ;

    R6a = k_POP.*POP .* f_CH4;
    R6b = k_POC.*POC .* f_CH4;
    R6c = k_DOP.*DOP .* f_CH4;
    R6d = k_DOC.*DOC .* f_CH4;
    R6f = k_Chl.*Chl .* f_CH4;



    Ra = R1a+R2a+R3a+R4a+R5a+R6a;
    Rb = R1b+R2b+R3b+R4b+R5b+R6b;
    Rc = R1c+R2c+R3c+R4c+R5c+R6c;
    Rd = R1d+R2d+R3d+R4d+R5d+R6d;
    Rf = R1f+R2f+R3f+R4f+R5f+R6f;

    R1 = R1a+R1b+R1c+R1d+R1f;
    R2 = R2a+R2b+R2c+R2d+R2f;
    R3 = R3a+R3b+R3c+R3d+R3f;
    R4 = R4a+R4b+R4c+R4d+R4f;
    R5 = R5a+R5b+R5c+R5d+R5f;
    R6 = R6a+R6b+R6c+R6d+R6f;


    R11 = k_tsox * Ox .* Sum_H2S;
    R12 = k_tS_Fe * FeOH3 .*  Sum_H2S;

    R13 = k_Feox .* Fe2 .* Ox;
    % NOTE: Due to the reaction is too fast and could cause overshooting:
    % we need to make this check if R*dt > Conc of source:
    % if R*dt > Conc then R13 = C/dt
    % if R*dt < Conc then R13 = R13
    % R13 = (R13.*dt < Fe2/50).*R13 + (R13.*dt > Fe2/50).* R13 ./ 1000;
    % R13 = (R13.*dt < Fe2).*R13 + (R13.*dt > Fe2).* Fe2 ./ (dt) * 0.5;
    % R13 = (R13.*dt < Ox).*R13 + (R13.*dt > Ox).* Ox ./ (dt) * 0.5;

    % R14 = k_amox * Ox ./ (Km_oxao + Ox) .* (NH4 ./ (Km_amao + NH4)); % NOTE: Doesnt work - Highly unstable.
    R14 = k_amox  .* NH4 .* Ox;
    % R14 = (R14.*dt < NH4).*R14 + (R14.*dt > NH4).* NH4 ./ (dt) * 0.5;
    % R14 = (R14.*dt < Ox).*R14 + (R14.*dt > Ox).* Ox ./ (dt) * 0.5;

    R21a = k_oms * Sum_H2S .* POP;
    R21b = k_oms * Sum_H2S .* POC;
    R21c = k_oms * Sum_H2S .* DOP;
    R21d = k_oms * Sum_H2S .* DOC;
    R21f = k_oms * Sum_H2S .* Chl;

    R22 = k_FeSpre .* FeS .* S0;
    R23 = k_rhom * Ox .* FeS;
    R24 = k_FeS2pre .* FeS .* Sum_H2S;

    % NOTE: Could cause instability. These rates are too high when pH > 7
    R25a = k_Fe_pre .* Fe2 .* HS .* (Sat_FeS > 1);
    R25b  =  k_Fe_dis .* FeS .* (Sat_FeS < 1);
    % R25a = (R25a >= 0) .* R25a; % can only be non negative
    % R25b  = (R25b >= 0) .* R25b; % can only be non negative

    R26a = k_Spre * S0;
    R26b = k_Sdis .* S8;

    R31a = k_pdesorb_a * FeOH3 .* PO4;
    R31b = 4 * (Cx2*R3a_P + Cx3*R3b_P + Cx2*R3c_P + Cx3*R3d_P + Cx1*R3f_P); % f_pfe .* (4 * R3 + 2 * R12);
    % R31b = (R31b.*dt < PO4adsa).*R31b + (R31b.*dt > PO4adsa).* PO4adsa ./ (dt) * 0.5;
    R32a = k_pdesorb_b * (FeOOH - PO4adsb) .* PO4;
    R32b = f_pfe .* (4 * R4);
    % R32b = (R32b.*dt < PO4adsb).*R32b + (R32b.*dt > PO4adsb).* PO4adsb ./ (dt) * 0.5;

    % R33 disabled now (no solid species Al=PO4)
    R33a = 0; % k_pdesorb_c .* PO4 .* AlOH3;
    R33b = 0;
    R34 = k_apa * (PO4 - kapa);
    R34 = (R34 >= 0) .* R34; % can only be non negative


    % saving rates
    if sediment_params.rate_estimator_switch
      r.R1a = R1a; r.R1b = R1b; r.R1c = R1c; r.R1d = R1d; ; r.R1f = R1f; r.R2a = R2a; r.R2b = R2b; r.R2c = R2c; r.R2d = R2d; ; r.R2f = R2f; r.R3a = R3a; r.R3b = R3b; r.R3c = R3c; r.R3d = R3d; ; r.R3f = R3f; r.R4a = R4a; r.R4b = R4b; r.R4c = R4c; r.R4d = R4d; ; r.R4f = R4f; r.R5a = R5a; r.R5b = R5b; r.R5c = R5c; r.R5d = R5d; r.R5f = R5f;  r.R6a = R6a; r.R6b = R6b; r.R6c = R6c; r.R6d = R6d; ; r.R6f = R6f; r.Ra = Ra; r.Rb = Rb; r.Rc = Rc; r.Rd = Rd; ; r.Rf = Rf; r.R1 = R1; r.R2 = R2; r.R3 = R3; r.R4 = R4; r.R5 = R5; r.R6 = R6;  r.R11 = R11; r.R12 = R12; r.R13  = R13; r.R14 = R14; r.R21a = R21a; r.R21b = R21b; r.R21c = R21c; r.R21d = R21d; ; r.R21f = R21f; r.R21f = R21f; r.R22 = R22; r.R23 = R23; r.R24 = R24; r.R25a = R25a; r.R25b  = R25b; r.R26a = R26a; r.R26b = R26b; r.R31a = R31a; r.R31b  = R31b; r.R32a = R32a; r.R32b = R32b; r.R33a = R33a; r.R33b = R33b; r.R34 = R34;
    else
      r = 0;
    end

    % for stoichiometry check:
    % Canavan, R. W., Slomp, C. P., Jourabchi, P., Van Cappellen, P., Laverman, A. M., & van den Berg, G. A. (2006). Organic matter mineralization in sediment of a coastal freshwater lake and response to salinization. Geochimica Et Cosmochimica Acta, 70(11), 2836–2855. http://doi.org/10.1016/j.gca.2006.03.012


    % F = 1./phi;
    F = (1-phi) ./ phi;

    dcdt(:,1)  = - bioirrigation(Ox, alfax, phi) +  -0.25 * R13  - 2 * R14  - (Cx2*R1a + Cx3*R1b+Cx1*R1f) .* F - (Cx2*R1c + Cx3*R1d) - 3 * R23; % Ox
    dcdt(:,2)  = -Ra - R21a; % POP
    dcdt(:,3)  = -Rb - R21b; % POC
    dcdt(:,4)  = - bioirrigation(NO3, alfax, phi) +  - 0.8*(Cx2*R2a+Cx2*R2b++Cx1*R2f) .* F - 0.8*(Cx2*R2c+Cx2*R2d)+ R14; % NO3
    dcdt(:,5)  = -4 * (Cx2*R3a_Fe + Cx3*R3b_Fe + Cx2*R3c_Fe + Cx3*R3d_Fe+ Cx1*R3f_Fe) - 2*R12 + R13./ F - R31a; % FeOH3
    dcdt(:,6)  = - bioirrigation(SO4, alfax, phi) +  - 0.5*(Cx2*R5a + Cx3*R5b+ Cx1*R5f) .* F -0.5*(Cx2*R5c + Cx3*R5d)+ R11; % SO4
    dcdt(:,7)  = - bioirrigation(NH4, alfax, phi) +  (Ny2 * Ra + Ny3 * Rb+ Ny1 * Rf) .* F + (Ny2 * Rc + Ny3 * Rd) - R14; % NH4
    dcdt(:,8)  = - bioirrigation(Fe2, alfax, phi) +  4*(Cx2*R3a + Cx3*R3b+ Cx1*R3f) .* F + 4* (Cx2*R3c + Cx3*R3d) + 4*(Cx2*R4a + Cx3*R4b+ Cx1*R4f) .* F + 4 * (Cx2*R4c + Cx3*R4d) + 2*R12 - R13 + R25b - R25a; % Fe2
    dcdt(:,9)  = -4*(Cx2*R4a + Cx3*R4b + Cx2*R4c + Cx3*R4d+ Cx1*R4f) + R23; % FeOOH
    dcdt(:,10) = - bioirrigation(H2S, alfax, phi); % H2S
    dcdt(:,11) = - bioirrigation(HS, alfax, phi) +  0.5*(Cx2*R5a + Cx3*R5b + Cx1*R5f) .* F + 0.5 * (Cx2*R5c + Cx3*R5d) - R11 - R12 + R25b - R25a - R21a - R21b - R21c - R21d - R21f -R24; % HS
    dcdt(:,12) =  - R22 - 4*R23 -R24 + R25a - R25b ; % FeS
    dcdt(:,13) = - R22 - R26a + R12 + R26b; % S0
    dcdt(:,14) = - bioirrigation(PO4, alfax, phi) +  (Pz2 * Ra + Pz3 * Rb + Pz1 * Rf) .* F + (Pz2 * Rc + Pz3 * Rd) + R31b + R32b - 2 * R34 - R33a - R31a - R32a; % PO4
    dcdt(:,15) = 4*R23 - R26b + R26a; % S8
    dcdt(:,16) = + R22 + R24; % FeS2
    dcdt(:,17) = -R33a; % AlOH3
    dcdt(:,18) = R31a - R31b; % PO4adsa
    dcdt(:,19) = R32a - R32b; % PO4adsb
    dcdt(:,20) = - bioirrigation(Ca2, alfax, phi) -3*R34; % Ca2
    dcdt(:,21) = R34; % Ca3PO42
    dcdt(:,22) = R21a + R21b + R21c + R21d + R21f; % OMS
    dcdt(:,23) = 0; % H
    dcdt(:,24) = 0; % OH
    dcdt(:,25) = - bioirrigation(CO2, alfax, phi)  +  ((Cx2 - Ny2 + 2*Pz2)*R1a + (Cx3 - Ny3 + 2*Pz3)*R1b + (Cx1 - Ny1 + 2*Pz1)*R1f + (0.2*Cx2 - Ny2 + 2*Pz2)*R2a +  (0.2*Cx3 - Ny3 + 2*Pz3)*R2b +  (0.2*Cx1 - Ny1 + 2*Pz1)*R2f - (7*Cx2 + Ny2 + 2*Pz2)*(R3a+R4a) - (7*Cx3 + Ny3 + 2*Pz3)*(R3b+R4b) - (7*Cx1 + Ny1 + 2*Pz1)*(R3f+R4f)  - (Ny2 - 2*Pz2)*R5a + (Ny3 - 2*Pz3)*R5b + (Ny1 - 2*Pz1)*R5f) .* F  +  (Cx2 - Ny2 + 2*Pz2)*R1c + (Cx3 - Ny3 + 2*Pz3)*R1d + (0.2*Cx2 - Ny2 + 2*Pz2)*R2c +  (0.2*Cx3 - Ny3 + 2*Pz3)*R2d - (7*Cx2 + Ny2 + 2*Pz2)*(R3c+R4c) - (7*Cx3 + Ny3 + 2*Pz3)*(R3d+R4d)  - (Ny2 - 2*Pz2)*R5c + (Ny3 - 2*Pz3)*R5d + 2*R13 + 2*R14;  % CO2
    dcdt(:,26) = - bioirrigation(CO3, alfax, phi) ; % CO3
    dcdt(:,27) = - bioirrigation(HCO3, alfax, phi) +  ((0.8*Cx2 + Ny2 - 2*Pz2)*R2a + (0.8*Cx3 + Ny3 - 2*Pz3)*R2b + (0.8*Cx1 + Ny1 - 2*Pz1)*R2f  + (8*Cx2+Ny2-2*Pz2)*(R3a + R4a) +(8*Cx3+Ny3-2*Pz3)*(R3b + R4b) +(8*Cx1+Ny1-2*Pz1)*(R3f + R4f)  + (Cx2+Ny2-2*Pz2)*R5a + (1*Cx3+Ny3-2*Pz3)*R5b + (1*Cx1+Ny1-2*Pz1)*R5f ) .* F + (0.8*Cx2 + Ny2 - 2*Pz2)*R2c + (0.8*Cx3 + Ny3 - 2*Pz3)*R2d + (8*Cx2+Ny2-2*Pz2)*(R3c + R4c) +(8*Cx3+Ny3-2*Pz3)*(R3d + R4d) + (Cx2+Ny2-2*Pz2)*R5c + (1*Cx3+Ny3-2*Pz3)*R5d -  2*R13 - 2*R14; % HCO3
    dcdt(:,28) = - bioirrigation(NH3, alfax, phi) ; % NH3
    dcdt(:,29) = - bioirrigation(H2CO3, alfax, phi) ; % H2CO3
    dcdt(:,30) = - bioirrigation(DOP, alfax, phi)  - Rc - R21c; % DOP
    dcdt(:,31) = - bioirrigation(DOC, alfax, phi)  - Rd - R21d; % DOC
    dcdt(:,32) = -Rf - R21f; ; % Chl
    dcdt(:,33) = - bioirrigation(CH4, alfax, phi) + F.*R6 ; %+ R16 - R17 ;  % CH4
end

