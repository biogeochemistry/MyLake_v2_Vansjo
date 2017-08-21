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
  CH4aq(:,1) = sediment_concentrations.CH4aq;
  CH4g(:,1) = sediment_concentrations.CH4g;


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
  D_CH4aq = sediment_params.D_CH4aq;
  D_CH4g = sediment_params.D_CH4g;
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
  [CH4aq_AL, CH4aq_AR] = sediment_matrix_templates{21,1:2};
  [CH4g_AL, CH4g_AR] = sediment_matrix_templates{22,1:2};


  % Solving equations!!!
  % =========================================================================================================

  for i=2:m
    % =======================================================================================================
    % Solving Reaction eq-s
    % =======================================================================================================
    C0 = [O2(:,i-1), POP(:,i-1), POC(:,i-1), NO3(:,i-1), FeOH3(:,i-1), SO4(:,i-1), NH4(:,i-1), Fe2(:,i-1), FeOOH(:,i-1), H2S(:,i-1), HS(:,i-1), FeS(:,i-1), S0(:,i-1), PO4(:,i-1), S8(:,i-1), FeS2(:,i-1), AlOH3(:,i-1), PO4adsa(:,i-1), PO4adsb(:,i-1), Ca2(:,i-1), Ca3PO42(:,i-1), OMS(:,i-1), H(:,i-1), OH(:,i-1), CO2(:,i-1), CO3(:,i-1), HCO3(:,i-1), NH3(:,i-1), H2CO3(:,i-1), DOP(:,i-1), DOC(:,i-1), Chl(:,i-1), CH4aq(:,i-1), CH4g(:,i-1)];

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
      CH4aq(:,i-1)   = C_new(:,33);
      CH4g(:,i-1)   = C_new(:,34);


    % =======================================================================================================
    % Solving Transport eq-s
    % =======================================================================================================

      O2(:,i) = pde_solver_dirichlet(O2_AL, O2_AR, O2(:,i-1), sediment_bc.O2_c);
      POP(:,i) = pde_solver_neumann(POP_AL, POP_AR, POP(:,i-1), sediment_bc.POP_fx, sediment_params.solid_flux_coef);
      POC(:,i) = pde_solver_neumann(POC_AL, POC_AR, POC(:,i-1), sediment_bc.POC_fx, sediment_params.solid_flux_coef);
      NO3(:,i) = pde_solver_dirichlet(NO3_AL, NO3_AR, NO3(:,i-1), sediment_bc.NO3_c);
      FeOH3(:,i) = pde_solver_neumann(FeOH3_AL, FeOH3_AR, FeOH3(:,i-1), sediment_bc.FeOH3_fx, sediment_params.solid_flux_coef);
      SO4(:,i) = pde_solver_dirichlet(SO4_AL, SO4_AR, SO4(:,i-1), sediment_bc.SO4_c);
      NH4(:,i) = pde_solver_dirichlet(NH4_AL, NH4_AR, NH4(:,i-1), sediment_bc.NH4_c);
      Fe2(:,i) = pde_solver_dirichlet(Fe2_AL, Fe2_AR, Fe2(:,i-1), sediment_bc.Fe2_c);
      FeOOH(:,i) = pde_solver_neumann(FeOOH_AL, FeOOH_AR, FeOOH(:,i-1), sediment_bc.FeOOH_fx, sediment_params.solid_flux_coef);
      H2S(:,i) = pde_solver_dirichlet(H2S_AL, H2S_AR, H2S(:,i-1), sediment_bc.H2S_c);
      HS(:,i) = pde_solver_dirichlet(HS_AL, HS_AR, HS(:,i-1), sediment_bc.HS_c);
      FeS(:,i) = pde_solver_neumann(FeS_AL, FeS_AR, FeS(:,i-1), sediment_bc.FeS_fx, sediment_params.solid_flux_coef);
      S0(:,i) = pde_solver_dirichlet(S0_AL, S0_AR, S0(:,i-1), sediment_bc.S0_c);
      PO4(:,i) = pde_solver_dirichlet(PO4_AL, PO4_AR, PO4(:,i-1), sediment_bc.PO4_c);
      S8(:,i) = pde_solver_neumann(S8_AL, S8_AR, S8(:,i-1), sediment_bc.S8_fx, sediment_params.solid_flux_coef);
      FeS2(:,i) = pde_solver_neumann(FeS2_AL, FeS2_AR, FeS2(:,i-1), sediment_bc.FeS2_fx, sediment_params.solid_flux_coef);
      AlOH3(:,i) = pde_solver_neumann(AlOH3_AL, AlOH3_AR, AlOH3(:,i-1), sediment_bc.AlOH3_fx, sediment_params.solid_flux_coef);
      PO4adsa(:,i) = pde_solver_neumann(PO4adsa_AL, PO4adsa_AR, PO4adsa(:,i-1), sediment_bc.PO4adsa_fx, sediment_params.solid_flux_coef);
      PO4adsb(:,i) = pde_solver_neumann(PO4adsb_AL, PO4adsb_AR, PO4adsb(:,i-1), sediment_bc.PO4adsb_fx, sediment_params.solid_flux_coef);
      Ca2(:,i) = pde_solver_dirichlet(Ca2_AL, Ca2_AR, Ca2(:,i-1), sediment_bc.Ca2_c);
      Ca3PO42(:,i) = pde_solver_neumann(Ca3PO42_AL, Ca3PO42_AR, Ca3PO42(:,i-1), sediment_bc.Ca3PO42_fx, sediment_params.solid_flux_coef);
      OMS(:,i) = pde_solver_neumann(OMS_AL, OMS_AR, OMS(:,i-1), sediment_bc.OMS_fx, sediment_params.solid_flux_coef);
      H(:,i) = pde_solver_dirichlet(H_AL, H_AR, H(:,i-1), sediment_bc.H_c);
      OH(:,i) = pde_solver_dirichlet(OH_AL, OH_AR, OH(:,i-1), sediment_bc.OH_c);
      CO2(:,i) = pde_solver_dirichlet(CO2_AL, CO2_AR, CO2(:,i-1), sediment_bc.CO2_c);
      CO3(:,i) = pde_solver_dirichlet(CO3_AL, CO3_AR, CO3(:,i-1), sediment_bc.CO3_c);
      HCO3(:,i) = pde_solver_dirichlet(HCO3_AL, HCO3_AR, HCO3(:,i-1), sediment_bc.HCO3_c);
      NH3(:,i) = pde_solver_dirichlet(NH3_AL, NH3_AR, NH3(:,i-1), sediment_bc.NH3_c);
      H2CO3(:,i) = pde_solver_dirichlet(H2CO3_AL, H2CO3_AR, H2CO3(:,i-1), sediment_bc.H2CO3_c);
      DOP(:,i) = pde_solver_dirichlet(DOP_AL, DOP_AR, DOP(:,i-1), sediment_bc.DOP_c);
      DOC(:,i) = pde_solver_dirichlet(DOC_AL, DOC_AR, DOC(:,i-1), sediment_bc.DOC_c);
      Chl(:,i) = pde_solver_neumann(Chl_AL, Chl_AR, Chl(:,i-1), sediment_bc.Chl_fx, sediment_params.solid_flux_coef);
      CH4aq(:,i) = pde_solver_dirichlet(CH4aq_AL, CH4aq_AR, CH4aq(:,i-1), sediment_bc.CH4aq_c);
      CH4g(:,i) = pde_solver_dirichlet(CH4g_AL, CH4g_AR, CH4g(:,i-1), sediment_bc.CH4g_fx);

      % Estimate fluxes:

      sediment_bioirrigation_fluxes.O2(i-1)   = integrate_over_depth_2( bioirrigation(O2(:, i), alfax, phi), x);
      sediment_bioirrigation_fluxes.PO4(i-1)  = integrate_over_depth_2( bioirrigation(PO4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.Fe2(i-1)  = integrate_over_depth_2( bioirrigation(Fe2(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.NO3(i-1)  = integrate_over_depth_2( bioirrigation(NO3(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.NH4(i-1)  = integrate_over_depth_2( bioirrigation(NH4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.SO4(i-1)  = integrate_over_depth_2( bioirrigation(SO4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.DOP(i-1) = integrate_over_depth_2( bioirrigation(DOP(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.DOC(i-1) = integrate_over_depth_2( bioirrigation(DOC(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CH4aq(i-1) = integrate_over_depth_2( bioirrigation(CH4aq(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CH4g(i-1) = integrate_over_depth_2( bioirrigation(CH4g(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CO2(i-1) = integrate_over_depth_2( bioirrigation(CO2(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.HCO3(i-1) = integrate_over_depth_2( bioirrigation(HCO3(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CO3(i-1) = integrate_over_depth_2( bioirrigation(CO3(:, i),  alfax,  phi), x);

      sediment_transport_fluxes.POP(i-1)          = -sediment_bc.POP_fx; % * (1-phi(1)) ./ phi(i);
      sediment_transport_fluxes.Chl(i-1)          = -sediment_bc.Chl_fx; % * (1-phi(1)) ./ phi(i);
      sediment_transport_fluxes.POC(i-1)          = -sediment_bc.POC_fx; % * (1-phi(1)) ./ phi(i);
      sediment_transport_fluxes.FeOH3(i-1)        = -sediment_bc.FeOH3_fx; % * (1-phi(1)) ./ phi(i);
      sediment_transport_fluxes.AlOH3(i-1)        = -sediment_bc.AlOH3_fx; % * (1-phi(1)) ./ phi(i);
      sediment_transport_fluxes.PO4adsa(i-1)      = -sediment_bc.PO4adsa_fx; % * (1-phi(1)) ./ phi(i);
      sediment_transport_fluxes.PO4adsb(i-1)      = -sediment_bc.PO4adsb_fx; % * (1-phi(1)) ./ phi(i);
      sediment_transport_fluxes.O2(i-1)           = top_sediment_diffusion_flux(O2(:, i), D_O2, dx, phi) + top_sediment_advection_flux(O2(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.PO4(i-1)          = top_sediment_diffusion_flux(PO4(:, i), D_PO4, dx, phi) + top_sediment_advection_flux(PO4(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.NO3(i-1)          = top_sediment_diffusion_flux(NO3(:, i), D_NO3, dx, phi) + top_sediment_advection_flux(NO3(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.Fe2(i-1)          = top_sediment_diffusion_flux(Fe2(:, i), D_Fe2, dx, phi) + top_sediment_advection_flux(Fe2(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.NH4(i-1)          = top_sediment_diffusion_flux(NH4(:, i), D_NH4, dx, phi) + top_sediment_advection_flux(NH4(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.SO4(i-1)          = top_sediment_diffusion_flux(SO4(:, i), D_SO4, dx, phi) + top_sediment_advection_flux(SO4(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.DOP(i-1)         = top_sediment_diffusion_flux(DOP(:, i), D_DOP, dx, phi) + top_sediment_advection_flux(DOP(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.DOC(i-1)         = top_sediment_diffusion_flux(DOC(:, i), D_DOC, dx, phi) + top_sediment_advection_flux(DOC(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.CH4aq(i-1)         = top_sediment_diffusion_flux(CH4aq(:, i), D_DOC, dx, phi) + top_sediment_advection_flux(CH4aq(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.CH4g(i-1)         = top_sediment_diffusion_flux(CH4g(:, i), D_DOC, dx, phi) + top_sediment_advection_flux(CH4g(:, i), sediment_params.w_CH4, phi);  % NOTE: rising velocity of methane
      sediment_transport_fluxes.CO2(i-1)         = top_sediment_diffusion_flux(CO2(:, i), D_DOC, dx, phi) + top_sediment_advection_flux(CO2(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.HCO3(i-1)         = top_sediment_diffusion_flux(HCO3(:, i), D_DOC, dx, phi) + top_sediment_advection_flux(HCO3(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.CO3(i-1)         = top_sediment_diffusion_flux(CO3(:, i), D_DOC, dx, phi) + top_sediment_advection_flux(CO3(:, i), sediment_params.w, phi);




    % pH Module
    if sediment_params.pH_algorithm ~= 0
      [H(:,i), OH(:,i), DOC(:,i), HCO3(:,i), CO2(:,i), CO3(:,i), NH3(:,i), NH4(:,i), HS(:,i), H2S(:,i)] = pH_module(sediment_params.pH_algorithm, H(:,i), OH(:,i), H2CO3(:,i), HCO3(:,i), CO2(:,i), CO3(:,i), NH3(:,i), NH4(:,i), HS(:,i), H2S(:,i), Fe2(:,i), Ca2(:,i), NO3(:,i), SO4(:,i), PO4(:,i), FeS(:,i), FeS2(:,i), FeOH3(:,i), FeOOH(:,i), Ca3PO42(:,i), PO4adsa(:,i), PO4adsb(:,i), sediment_bc.T, sediment_params);
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
  sediment_transport_fluxes.CH4aq          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CH4aq), 16042.5);
  sediment_transport_fluxes.CH4g          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CH4g), 16042.5);
  sediment_transport_fluxes.CO2          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CO2), 44009.5);
  sediment_transport_fluxes.HCO3          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.HCO3), 61016.8);
  sediment_transport_fluxes.CO3          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CO3), 60008.9);

  sediment_bioirrigation_fluxes.O2   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.O2), 31998);
  sediment_bioirrigation_fluxes.PO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.PO4), 30973.762);
  sediment_bioirrigation_fluxes.Fe2  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.Fe2), 55845);
  sediment_bioirrigation_fluxes.NO3  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.NO3), 62004);
  sediment_bioirrigation_fluxes.NH4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.NH4), 18038);
  sediment_bioirrigation_fluxes.SO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.SO4), 96062);
  sediment_bioirrigation_fluxes.DOP = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.DOP), 30973.762);
  sediment_bioirrigation_fluxes.DOC = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.DOC), 12010.7);
  sediment_bioirrigation_fluxes.CH4aq = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CH4aq), 16042.5);
  sediment_bioirrigation_fluxes.CH4g = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CH4g), 16042.5);
  sediment_bioirrigation_fluxes.CO2 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CO2), 44009.5);
  sediment_bioirrigation_fluxes.HCO3 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.HCO3), 61016.8);
  sediment_bioirrigation_fluxes.CO3 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CO3), 60008.9);


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
  sediment_concentrations.CH4aq = CH4aq(:,end);
  sediment_concentrations.CH4g = CH4g(:,end);

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

    if any(isnan(sediment_transport_fluxes.O2))| any(isnan(sediment_bc.POP_fx))| any(isnan(sediment_bc.POC_fx))| any(isnan(sediment_bc.FeOH3_fx))| any(isnan(O2)) | any(isnan(POP)) | any(isnan(POC)) | any(isnan(NO3)) | any(isnan(FeOH3)) | any(isnan(SO4)) | any(isnan(NH4)) | any(isnan(Fe2)) | any(isnan(FeOOH)) | any(isnan(H2S)) | any(isnan(HS)) | any(isnan(FeS)) | any(isnan(S0)) | any(isnan(PO4)) | any(isnan(S8)) | any(isnan(FeS2)) | any(isnan(AlOH3)) | any(isnan(PO4adsa)) | any(isnan(PO4adsb)) | any(isnan(H)) | any(isnan(Ca2)) | any(isnan(Ca3PO42)) | any(isnan(OMS)) | any(isnan(OH)) | any(isnan(HCO3)) | any(isnan(CO2)) | any(isnan(CO3)) | any(isnan(NH3)) | any(isnan(H2CO3)) | any(isnan(Chl)) | any(isnan(CH4aq)) | any(isnan(CH4g))
      error('Breaking out of Sediments function: NaN values');
    end


end



function [H, OH, H2CO3, HCO3, CO2, CO3, NH3, NH4, HS, H2S] = pH_module(algorithm, H, OH, H2CO3, HCO3, CO2, CO3, NH3, NH4, HS, H2S, Fe2, Ca2, NO3, SO4, PO4, FeS, FeS2, FeOH3, FeOOH, Ca3PO42, PO4adsa, PO4adsb, Temperature, sediment_params)
  %% pH_module: pH equilibrium function
  % 0. No pH module
  % 1. Phreeqc
  % 2. New algorithm by Markelov (under test)
  % 3. Phreeqc Py

    Kc1=5.01*10^(-7); Kc2=4.78*10^(-11); Knh=5.62*10^(-10); Khs=1.3*10^(-7); Kw=10^(-14); Kc0 = 1.7*10^(-3);

    if algorithm == 1 %
        in =[H HCO3 CO2 CO3 NH3 NH4 HS H2S OH H2CO3 Fe2 Ca2 NO3 SO4 PO4 FeS FeS2 FeOH3 FeOOH Ca3PO42 PO4adsa PO4adsb];
        [pH_est] = pH_phreeqc(size(H,1),in);


        H = 10.^(-pH_est');
        Ct = H2CO3 + HCO3 + CO3;
        a = alpha(-log10(H), [6.52, 10.56]);
        CO2 = Ct.*a(:,1);
        HCO3 = Ct.*a(:,2);
        CO3 = Ct.*a(:,3);
        Nt = NH3 + NH4;
        St = HS + H2S;
        OH = Kw./H;
        NH4 = (H./ (Knh + H)) .* Nt;
        NH3 = Nt - NH4;
        H2S = (H./ (Khs + H)) .* St;
        HS = St - H2S;
        H = H*1e3;
        OH = OH*1e3;

    elseif algorithm == 2
      for i=1:size(H,1)
        sediment_params.aq_system.carb_acid.conc = 1e-3*(CO2(i)+HCO3(i)+CO3(i));
        sediment_params.aq_system.amonia.conc = 1e-3*(NH4(i)+NH3(i));
        sediment_params.aq_system.sulf.conc = 1e-3*(H2S(i)+HS(i));
        sediment_params.aq_system.ca.conc = 1e-3*(Ca2(i));
        sediment_params.aq_system.fe2.conc = 1e-3*(Fe2(i));
        sediment_params.aq_system.no3.conc = 1e-3*(NO3(i));
        sediment_params.aq_system.so4.conc = 1e-3*(SO4(i));
        sediment_params.aq_system.p_acid.conc = 1e-3*(PO4(i));

        if i == 1
          pHs = linspace(0,14,1400)';
        else
          pHs = linspace(pHz(i-1)-0.5,pHz(i-1)+0.5,100)';
        end
        pHz(i) = new_pH_module(sediment_params.aq_system, pHs);

        res = bsxfun(@times, sediment_params.aq_system.carb_acid.conc, alpha(pHz(i), sediment_params.aq_system.carb_acid.pKs));
        CO2(i) = 1e3*res(1);
        HCO3(i) = 1e3*res(2);
        CO3(i) = 1e3*res(3);
        res = bsxfun(@times, sediment_params.aq_system.amonia.conc, alpha(pHz(i), sediment_params.aq_system.amonia.pKs));
        NH4(i) = 1e3*res(1);
        NH3(i) = 1e3*res(2);
        res = bsxfun(@times, sediment_params.aq_system.sulf.conc, alpha(pHz(i), sediment_params.aq_system.sulf.pKs));
        H2S(i) = 1e3*res(1);
        HS(i) = 1e3*res(2);
        H(i) = 10^(-pHz(i))*1e3;
        OH(i) = 10^(-14+pHz(i))*1e3;
      end
    end
end


function C_new = pde_solver_dirichlet(AL, AR, C_old, const_bc)
    C_old(1) = const_bc;
    temp = AR*C_old;
    % temp(1) = const_bc;
    C_new = AL\ temp;
    % C_new(1) = const_bc;
    C_new = (C_new>0).*C_new;
end


function C_new = pde_solver_neumann(AL, AR, C_old, flux_bc, coef)
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

function [flux] = top_sediment_advection_flux(C, w, phi)
  % calculates flux of the particular dissolved specie through the top boundary of the sediment
  % in [ mg m-2 d-1 ] units
  % C - concentration
  % phi - porosity (no porosity because C is the concentration in pores (not bulk))
  % w - advection velocity
  % minus sign because of velocity signs

  % fourth-order
  flux = - phi(1) * w * C(1) ;  %  [umol/cm^2/y]
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

    O2      = C(:,1) .* (C(:,1)>0) ;
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
    CH4aq    = C(:,33) .* (C(:,33)>0) ;
    CH4g    = C(:,34) .* (C(:,34)>0) ;

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
    k_fesox   = sediment_params.k_fesox;
    k_tS_Fe = sediment_params.k_tS_Fe;
    Ks_FeS = sediment_params.Ks_FeS;
    k_Fe_dis = sediment_params.k_Fe_dis;
    k_Fe_pre = sediment_params.k_Fe_pre;
    k_apa  = sediment_params.k_apa;
    kapa = sediment_params.kapa;
    k_oms = sediment_params.k_oms;
    k_tsox = sediment_params.k_tsox;
    k_FeSpre = sediment_params.k_FeSpre;
    k_ch4_o2 = sediment_params.k_ch4_o2;
    k_ch4_so4 = sediment_params.k_ch4_so4;
    CH4_solubility = sediment_params.CH4_solubility;
    CO2_solubility = sediment_params.CO2_solubility;
    k_ch4_dis = sediment_params.k_ch4_dis;
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
    tot_FeOOH = PO4adsb + FeOOH;

    f_O2    = O2 ./  (Km_O2 + O2);
    f_NO3   = NO3 ./  (Km_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_FeOH3 = tot_FeOH3 ./  (Km_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_FeOOH = tot_FeOOH ./  (Km_FeOOH + tot_FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_SO4   = SO4 ./ (Km_SO4 + SO4 ) .* Kin_FeOOH ./ (Kin_FeOOH + FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_CH4 = 1 - f_O2 - f_NO3 - f_FeOH3 - f_FeOOH - f_SO4;
    f_CH4 = f_CH4.*(f_CH4>0);

    % % Solubility of CH4aq and switch function as in Canavan
    % delta = (CH4aq - CH4_solubility/2)/Sm;
    % fm = (tanh(-delta) + 1)/2;








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


    part_PO4ads_tot_FeOH3 = PO4adsa ./ (tot_FeOH3+1e-16); % ratio of ads P to total Fe(III)
    R3a_Fe = (1 - part_PO4ads_tot_FeOH3) .* k_POP.*POP .* f_FeOH3;
    R3b_Fe = (1 - part_PO4ads_tot_FeOH3) .* k_POC .*POC .* f_FeOH3;
    R3c_Fe = (1 - part_PO4ads_tot_FeOH3) .* k_DOP .*DOP .* f_FeOH3;
    R3d_Fe = (1 - part_PO4ads_tot_FeOH3) .* k_DOC .*DOC .* f_FeOH3;
    R3f_Fe = (1 - part_PO4ads_tot_FeOH3) .* k_Chl .*Chl .* f_FeOH3;
    R3a_P = part_PO4ads_tot_FeOH3 .* k_POP.*POP .* f_FeOH3;
    R3b_P = part_PO4ads_tot_FeOH3 .* k_POC .*POC .* f_FeOH3;
    R3c_P = part_PO4ads_tot_FeOH3 .* k_DOP .*DOP .* f_FeOH3;
    R3d_P = part_PO4ads_tot_FeOH3 .* k_DOC .*DOC .* f_FeOH3;
    R3f_P = part_PO4ads_tot_FeOH3 .* k_Chl .*Chl .* f_FeOH3;
    R3a = R3a_Fe + R3a_P;
    R3b = R3b_Fe + R3b_P;
    R3c = R3c_Fe + R3c_P;
    R3d = R3d_Fe + R3d_P;
    R3f = R3f_Fe + R3f_P;


    part_PO4ads_tot_FeOOH = PO4adsb ./ (tot_FeOOH+1e-16); % ratio of ads P to total Fe(III)
    R4a_Fe = (1 - part_PO4ads_tot_FeOOH) .* k_POP.*POP .* f_FeOOH;
    R4b_Fe = (1 - part_PO4ads_tot_FeOOH) .* k_POC .*POC .* f_FeOOH;
    R4c_Fe = (1 - part_PO4ads_tot_FeOOH) .* k_DOP .*DOP .* f_FeOOH;
    R4d_Fe = (1 - part_PO4ads_tot_FeOOH) .* k_DOC .*DOC .* f_FeOOH;
    R4f_Fe = (1 - part_PO4ads_tot_FeOOH) .* k_Chl .*Chl .* f_FeOOH;
    R4a_P = part_PO4ads_tot_FeOOH .* k_POP.*POP .* f_FeOOH;
    R4b_P = part_PO4ads_tot_FeOOH .* k_POC .*POC .* f_FeOOH;
    R4c_P = part_PO4ads_tot_FeOOH .* k_DOP .*DOP .* f_FeOOH;
    R4d_P = part_PO4ads_tot_FeOOH .* k_DOC .*DOC .* f_FeOOH;
    R4f_P = part_PO4ads_tot_FeOOH .* k_Chl .*Chl .* f_FeOOH;
    R4a = R4a_Fe + R4a_P;
    R4b = R4b_Fe + R4b_P;
    R4c = R4c_Fe + R4c_P;
    R4d = R4d_Fe + R4d_P;
    R4f = R4f_Fe + R4f_P;

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


    Sum_H2S = H2S + HS;;
    R11 = k_tsox * O2 .* Sum_H2S;
    R12 = k_tS_Fe * FeOH3 .*  Sum_H2S;

    R13 = k_Feox .* Fe2 .* O2;
    % NOTE: Due to the reaction is too fast and could cause overshooting:
    % we need to make this check if R*dt > Conc of source:
    % if R*dt > Conc then R13 = C/dt
    % if R*dt < Conc then R13 = R13
    % R13 = (R13.*dt < Fe2/50).*R13 + (R13.*dt > Fe2/50).* R13 ./ 1000;
    % R13 = (R13.*dt < Fe2).*R13 + (R13.*dt > Fe2).* Fe2 ./ (dt) * 0.5;
    % R13 = (R13.*dt < O2).*R13 + (R13.*dt > O2).* O2 ./ (dt) * 0.5;

    % R14 = k_amox * O2 ./ (Km_oxao + O2) .* (NH4 ./ (Km_amao + NH4)); % NOTE: Doesnt work - Highly unstable.
    R14 = k_amox  .* NH4 .* O2;
    % R14 = (R14.*dt < NH4).*R14 + (R14.*dt > NH4).* NH4 ./ (dt) * 0.5;
    % R14 = (R14.*dt < O2).*R14 + (R14.*dt > O2).* O2 ./ (dt) * 0.5;

    R15 = k_ch4_o2 .* CH4aq .* O2;
    R16 = k_ch4_so4 .* CH4aq .* SO4;


    CH4_over_sat = CH4aq - CH4_solubility;
    R17 = k_ch4_dis .* CH4_over_sat .* (CH4_over_sat > 0);

    CO2_over_sat = CO2 - CO2_solubility;


    R21a = k_oms * Sum_H2S .* POP;
    R21b = k_oms * Sum_H2S .* POC;
    R21c = k_oms * Sum_H2S .* DOP;
    R21d = k_oms * Sum_H2S .* DOC;
    R21f = k_oms * Sum_H2S .* Chl;

    R22 = k_FeSpre .* FeS .* S0;
    R23 = k_fesox * O2 .* FeS;
    R24 = k_FeS2pre .* FeS .* Sum_H2S;

    % NOTE: Could cause instability. These rates are too high when pH > 7
    Sat_FeS = Fe2*1e-3 .* Sum_H2S*1e-3 ./ (H*1e-3+1e-16).^2 ./ Ks_FeS;
    R25a = k_Fe_pre .* Fe2 .* Sum_H2S .* (Sat_FeS > 1); %(Sat_FeS-1)
    R25b  =  k_Fe_dis .* FeS .* (Sat_FeS < 1); % .* (1-Sat_FeS)
    % R25a = (R25a >= 0) .* R25a; % can only be non negative
    % R25b  = (R25b >= 0) .* R25b; % can only be non negative

    R26a = k_Spre * S0;
    R26b = k_Sdis .* S8;

    R31a = k_pdesorb_a * FeOH3 .* PO4;
    R31b = 4 * (Cx2*R3a_P + Cx3*R3b_P + Cx2*R3c_P + Cx3*R3d_P + Cx1*R3f_P); % f_pfe .* (4 * R3 + 2 * R12);
    % R31b = (R31b.*dt < PO4adsa).*R31b + (R31b.*dt > PO4adsa).* PO4adsa ./ (dt) * 0.5;
    R32a = k_pdesorb_b * FeOOH .* PO4;
    R32b = 4 * (Cx2*R4a_P + Cx3*R4b_P + Cx2*R4c_P + Cx3*R4d_P + Cx1*R4f_P);
    % R32b = (R32b.*dt < PO4adsb).*R32b + (R32b.*dt > PO4adsb).* PO4adsb ./ (dt) * 0.5;

    % R33 disabled now (no solid species Al=PO4)
    R33a = 0; % k_pdesorb_c .* PO4 .* AlOH3;
    R33b = 0;
    R34 = k_apa * (PO4 - kapa);
    R34 = (R34 >= 0) .* R34; % can only be non negative


    % saving rates
    if sediment_params.rate_estimator_switch
      r.R1a = R1a; r.R1b = R1b; r.R1c = R1c; r.R1d = R1d; ; r.R1f = R1f; r.R2a = R2a; r.R2b = R2b; r.R2c = R2c; r.R2d = R2d; ; r.R2f = R2f; r.R3a = R3a; r.R3b = R3b; r.R3c = R3c; r.R3d = R3d; ; r.R3f = R3f; r.R4a = R4a; r.R4b = R4b; r.R4c = R4c; r.R4d = R4d; ; r.R4f = R4f; r.R5a = R5a; r.R5b = R5b; r.R5c = R5c; r.R5d = R5d; r.R5f = R5f;  r.R6a = R6a; r.R6b = R6b; r.R6c = R6c; r.R6d = R6d; ; r.R6f = R6f; r.Ra = Ra; r.Rb = Rb; r.Rc = Rc; r.Rd = Rd; ; r.Rf = Rf; r.R1 = R1; r.R2 = R2; r.R3 = R3; r.R4 = R4; r.R5 = R5; r.R6 = R6;  r.R11 = R11; r.R12 = R12; r.R13  = R13; r.R14 = R14; r.R15 = R15; r.R16 = R16; r.R21a = R21a; r.R21b = R21b; r.R21c = R21c; r.R21d = R21d; ; r.R21f = R21f; r.R21f = R21f; r.R22 = R22; r.R23 = R23; r.R24 = R24; r.R25a = R25a; r.R25b  = R25b; r.R26a = R26a; r.R26b = R26b; r.R31a = R31a; r.R31b  = R31b; r.R32a = R32a; r.R32b = R32b; r.R33a = R33a; r.R33b = R33b; r.R34 = R34;
    else
      r = 0;
    end

    % for stoichiometry check:
    % Canavan, R. W., Slomp, C. P., Jourabchi, P., Van Cappellen, P., Laverman, A. M., & van den Berg, G. A. (2006). Organic matter mineralization in sediment of a coastal freshwater lake and response to salinization. Geochimica Et Cosmochimica Acta, 70(11), 2836–2855. http://doi.org/10.1016/j.gca.2006.03.012


    % F = 1./phi;
    F = (1-phi) ./ phi;

    dcdt(:,1)  = - bioirrigation(O2, alfax, phi) +  -0.25 * R13  - R15 - 2 * R14  - (Cx2*R1a + Cx3*R1b+Cx1*R1f) .* F - (Cx2*R1c + Cx3*R1d) - 3 * R23.*F ; % O2 (aq)
    dcdt(:,2)  = -Ra - R21a; % POP (solid)
    dcdt(:,3)  = -Rb - R21b; % POC (solid)
    dcdt(:,4)  = - bioirrigation(NO3, alfax, phi) +  - 0.8*(Cx2*R2a+Cx2*R2b+Cx1*R2f) .* F - 0.8*(Cx2*R2c+Cx2*R2d)+ R14; % NO3(aq)
    dcdt(:,5)  = -4 * (Cx2*R3a_Fe + Cx3*R3b_Fe + Cx2*R3c_Fe + Cx3*R3d_Fe+ Cx1*R3f_Fe) - 2*R12 + R13./ F - R31a; % FeOH3(solid)
    dcdt(:,6)  = - bioirrigation(SO4, alfax, phi) +  - 0.5*(Cx2*R5a + Cx3*R5b+ Cx1*R5f) .* F -0.5*(Cx2*R5c + Cx3*R5d)+ R11 - R16; % SO4(aq)
    dcdt(:,7)  = - bioirrigation(NH4, alfax, phi) +  (Ny2 * Ra + Ny3 * Rb+ Ny1 * Rf) .* F + (Ny2 * Rc + Ny3 * Rd) - R14; % NH4(aq)
    dcdt(:,8)  = - bioirrigation(Fe2, alfax, phi) +  4*(Cx2*R3a + Cx3*R3b+ Cx1*R3f) .* F + 4* (Cx2*R3c + Cx3*R3d) + 4*(Cx2*R4a + Cx3*R4b+ Cx1*R4f) .* F + 4 * (Cx2*R4c + Cx3*R4d) + 2*R12.*F - R13 + R25b.*F - R25a; % Fe2(aq)
    dcdt(:,9)  = -4*(Cx2*R4a_Fe + Cx3*R4b_Fe + Cx2*R4c_Fe + Cx3*R4d_Fe+ Cx1*R4f_Fe) + R23 - R32a; % FeOOH(solid)
    dcdt(:,10) = - bioirrigation(H2S, alfax, phi); % H2S(aq)
    dcdt(:,11) = - bioirrigation(HS, alfax, phi) +  0.5*(Cx2*R5a + Cx3*R5b + Cx1*R5f) .* F + 0.5 * (Cx2*R5c + Cx3*R5d) - R11 - R12.*F + R25b.*F - R25a + (R21a + R21b + R21c + R21d + R21f).*F -R24.*F + R16; % HS(aq)
    dcdt(:,12) =  - R22 - 4*R23 -R24 + R25a./F - R25b ; % FeS(solid)
    dcdt(:,13) = - R22.*F - R26a + R12.*F + R26b.*F; % S0(aq)
    dcdt(:,14) = - bioirrigation(PO4, alfax, phi) +  (Pz2 * Ra + Pz3 * Rb + Pz1 * Rf) .* F + (Pz2 * Rc + Pz3 * Rd) + R31b + R32b.*F + R32b.*F - R31a.*F - R32a.*F - R33a.*F - 2 * R34; % PO4(aq)
    dcdt(:,15) = 4*R23 - R26b + R26a./F; % S8(solid)
    dcdt(:,16) = + R22 + R24; % FeS2(solid)
    dcdt(:,17) = -R33a; % AlOH3(s)
    dcdt(:,18) = R31a - R31b; % PO4adsa(s)
    dcdt(:,19) = R32a - R32b; % PO4adsb(s)
    dcdt(:,20) = - bioirrigation(Ca2, alfax, phi) -3*R34; % Ca2(aq)
    dcdt(:,21) = R34./F; % Ca3PO42(s)
    dcdt(:,22) = R21a + R21b + R21c + R21d + R21f; % OMS(s)
    dcdt(:,23) = 0; % H(aq)
    dcdt(:,24) = 0; % OH(aq)
    dcdt(:,25) = - bioirrigation(CO2, alfax, phi) + (R15 +  ((Cx2 - Ny2 + 2*Pz2)*R1a + (Cx3 - Ny3 + 2*Pz3)*R1b + (Cx1 - Ny1 + 2*Pz1)*R1f + (0.2*Cx2 - Ny2 + 2*Pz2)*R2a +  (0.2*Cx3 - Ny3 + 2*Pz3)*R2b +  (0.2*Cx1 - Ny1 + 2*Pz1)*R2f  + (0.5 * Cx2 - Ny2 - 2*Pz2)*R6a + (0.5 * Cx3 - Ny3 - 2*Pz3)*R6b + (0.5 * Cx1 - Ny1 - 2*Pz1)*R6f) .* F  +  (Cx2 - Ny2 + 2*Pz2)*R1c + (Cx3 - Ny3 + 2*Pz3)*R1d + (0.2*Cx2 - Ny2 + 2*Pz2)*R2c +  (0.2*Cx3 - Ny3 + 2*Pz3)*R2d - (7*Cx2 + Ny2 + 2*Pz2)*(R3c+R4c) + (0.5*Cx2 - Ny2 + 2*Pz2)*R6c + (0.5*Cx3 - Ny3 + 2*Pz3)*R6d).*(CO2_over_sat<0) + (- (7*Cx2 + Ny2 + 2*Pz2)*(R3a+R4a) - (7*Cx3 + Ny3 + 2*Pz3)*(R3b+R4b) - (7*Cx1 + Ny1 + 2*Pz1)*(R3f+R4f)  - (Ny2 - 2*Pz2)*R5a - (Ny3 - 2*Pz3)*R5b - (Ny1 - 2*Pz1)*R5f) .* F  - (7*Cx3 + Ny3 + 2*Pz3)*(R3d+R4d)  - (Ny2 - 2*Pz2)*R5c - (Ny3 - 2*Pz3)*R5d - R16 + 2*R13 + 2*R14;  % CO2 (aq) NOTE: we need to add gas (Rename H2CO3 to gas)
    dcdt(:,26) = - bioirrigation(CO3, alfax, phi) ; % CO3(aq)
    dcdt(:,27) = - bioirrigation(HCO3, alfax, phi)+  ( - (Ny2 - 2*Pz2)*R1a - (Ny3 - 2*Pz3)*R1b - (Ny1 - 2*Pz1)*R1f +  (0.8*Cx2 + Ny2 - 2*Pz2)*R2a + (0.8*Cx3 + Ny3 - 2*Pz3)*R2b + (0.8*Cx1 + Ny1 - 2*Pz1)*R2f  + (8*Cx2+Ny2-2*Pz2)*(R3a + R4a) +(8*Cx3+Ny3-2*Pz3)*(R3b + R4b) +(8*Cx1+Ny1-2*Pz1)*(R3f + R4f)  + (Cx2+Ny2-2*Pz2)*R5a + (1*Cx3+Ny3-2*Pz3)*R5b + (1*Cx1+Ny1-2*Pz1)*R5f + (Ny2-2*Pz2)*R6a + (Ny3-2*Pz3)*R6b + (Ny1-2*Pz1)*R6f) .* F - (Ny2 - 2*Pz2)*R1c - (Ny3 - 2*Pz3)*R1d + (0.8*Cx2 + Ny2 - 2*Pz2)*R2c + (0.8*Cx3 + Ny3 - 2*Pz3)*R2d + (8*Cx2+Ny2-2*Pz2)*(R3c + R4c) +(8*Cx3+Ny3-2*Pz3)*(R3d + R4d) + (Cx2+Ny2-2*Pz2)*R5c + (1*Cx3+Ny3-2*Pz3)*R5d + (Ny2-2*Pz2)*R6c + (Ny3-2*Pz3)*R6d  -  2*R13 - 2*R14 + 2*R16; % HCO3(aq)
    dcdt(:,28) = - bioirrigation(NH3, alfax, phi) ; % NH3(aq)
    dcdt(:,29) = - bioirrigation(H2CO3, alfax, phi) ; % H2CO3 (aq)
    dcdt(:,30) = - bioirrigation(DOP, alfax, phi)  - Rc - R21c; % DOP (aq)
    dcdt(:,31) = - bioirrigation(DOC, alfax, phi)  - Rd - R21d; % DOC (aq)
    dcdt(:,32) = -Rf - R21f; ; % Chl (s)
    dcdt(:,33) = - bioirrigation(CH4aq, alfax, phi) + 0.5*(Cx2*R6a+Cx2*R6b+Cx1*R6f).* F .* (CH4_over_sat < 0) + 0.5*(Cx2*R6c + Cx3*R6d).* (CH4_over_sat < 0) - R15 - R16 - R17;  % CH4aq
    dcdt(:,34) = - bioirrigation(CH4g, alfax, phi) + 0.5*(Cx2*R6a+Cx2*R6b+Cx1*R6f).* F .* (CH4_over_sat > 0) + 0.5*(Cx2*R6c + Cx3*R6d).* (CH4_over_sat > 0) + R17;  % CH4g
end

