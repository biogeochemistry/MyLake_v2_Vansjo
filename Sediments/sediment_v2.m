function [ sediment_bioirrigation_fluxes, sediment_transport_fluxes, sediment_concentrations, sediment_additional_results] = sediment_v2(sediment_concentrations, sediment_params, sediment_matrix_templates, sediment_bc)
  % SEDIMENTS This function models the chemical process in the sediment

  global k_OM k_OMb k_DOM1 k_DOM2 Km_O2 Km_NO3 Km_FeOH3 Km_FeOOH Km_SO4 Km_oxao Km_amao Kin_O2 Kin_NO3  Kin_FeOH3 Kin_FeOOH k_amox k_Feox k_Sdis k_Spre k_FeS2pre k_pdesorb_c k_pdesorb_a k_pdesorb_b k_alum k_rhom   k_tS_Fe Ks_FeS k_Fe_dis k_Fe_pre k_apa  kapa k_oms k_tsox k_FeSpre f_pfe accel Cx1 Ny1 Pz1 Cx2 Ny2 Pz2 F Ny1 Ny2 Pz1 Pz2 alfax fi n


  Ox_prev = sediment_concentrations.Oxygen;
  OM_prev = sediment_concentrations.OM1;
  OMb_prev = sediment_concentrations.OM2;
  NO3_prev = sediment_concentrations.NO3;
  FeOH3_prev = sediment_concentrations.FeOH3;
  SO4_prev = sediment_concentrations.SO4;
  NH4_prev = sediment_concentrations.NH4;
  Fe2_prev = sediment_concentrations.Fe2;
  FeOOH_prev = sediment_concentrations.FeOOH;
  H2S_prev = sediment_concentrations.H2S;
  HS_prev  = sediment_concentrations.HS;
  FeS_prev = sediment_concentrations.FeS;
  S0_prev  = sediment_concentrations.S0;
  PO4_prev = sediment_concentrations.PO4;
  S8_prev = sediment_concentrations.S8;
  FeS2_prev = sediment_concentrations.FeS2;
  AlOH3_prev = sediment_concentrations.AlOH3;
  PO4adsa_prev = sediment_concentrations.PO4adsa;
  PO4adsb_prev = sediment_concentrations.PO4adsb;
  Ca2_prev = sediment_concentrations.Ca2;
  Ca3PO42_prev = sediment_concentrations.Ca3PO42;
  OMS_prev = sediment_concentrations.OMS;
  OH_prev = sediment_concentrations.OH;
  CO2_prev = sediment_concentrations.CO2;
  CO3_prev = sediment_concentrations.CO3;
  HCO3_prev = sediment_concentrations.HCO3;
  NH3_prev = sediment_concentrations.NH3;
  H_prev = sediment_concentrations.H;
  H2CO3_prev = sediment_concentrations.H2CO3;
  DOM1_prev = sediment_concentrations.DOM1;
  DOM2_prev = sediment_concentrations.DOM2;


  % model domain:
  n  = sediment_params.n; %points in spatial grid
  depth = sediment_params.depth; %sediment depth
  years = sediment_params.years; %1 day
  ts    = sediment_params.ts; % time step
  Db    = sediment_params.Db; % is effective diffusion due to bioturbation, Canavan et al D_bio between 0-5, 5 in the top layers
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
  D_DOM1= sediment_params.D_DOM1;
  D_DOM2= sediment_params.D_DOM2;
  Cx1   = sediment_params.Cx1;
  Ny1   = sediment_params.Ny1;
  Pz1   = sediment_params.Pz1;
  Cx2   = sediment_params.Cx2;
  Ny2   = sediment_params.Ny2;
  Pz2   = sediment_params.Pz2;
  fi    = sediment_params.fi;
  % TODO:we need to reestimate F due to non-constant porosity profile. For this we need rhob = solid phase density
  F     = sediment_params.F;  %conversion factor = rhob * (1-fi) / fi ; where fi = porosity and rhob = solid phase density




  k_OM =  sediment_params.k_OM;
  k_OMb = sediment_params.k_OMb;
  k_DOM1 = sediment_params.k_DOM1;
  k_DOM2 = sediment_params.k_DOM2;
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
  [OM_AL, OM_AR]= sediment_matrix_templates{1,1:2};
  [OMb_AL, OMb_AR]= sediment_matrix_templates{1,1:2};
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

  % Solute species:
  [Ox_AL, Ox_AR] = sediment_matrix_templates{2,1:2};
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
  [DOM1_AL, DOM1_AR] = sediment_matrix_templates{19,1:2};
  [DOM2_AL, DOM2_AR] = sediment_matrix_templates{20,1:2};



  % Allocation of the memory for concentration with initial condition: (umol/cm3(aq)) or (umol/cm3(solid))
  % ===============================================================
  % NOTE: All columns consist of init vector of concentrations

  Ox = zeros(n,m);
  Ox(:,1) = Ox_prev;
  OM = zeros(n,m);
  OM(:,1) = OM_prev;
  OMb = zeros(n,m);
  OMb(:,1) = OMb_prev;
  NO3 = zeros(n,m);
  NO3(:,1) = NO3_prev;
  FeOH3 = zeros(n,m);
  FeOH3(:,1) = FeOH3_prev;
  SO4 = zeros(n,m);
  SO4(:,1) = SO4_prev;
  Fe2 = zeros(n,m);
  Fe2(:,1) = Fe2_prev;
  FeOOH = zeros(n,m);
  FeOOH(:,1) = FeOOH_prev;
  FeS = zeros(n,m);
  FeS(:,1) = FeS_prev;
  S0 = zeros(n,m);
  S0(:,1) = S0_prev;
  PO4 = zeros(n,m);
  PO4(:,1) = PO4_prev;
  S8 = zeros(n,m);
  S8(:,1) = S8_prev;
  FeS2 = zeros(n,m);
  FeS2(:,1) = FeS2_prev;
  AlOH3 = zeros(n,m);
  AlOH3(:,1) = AlOH3_prev;
  PO4adsa = zeros(n,m);
  PO4adsa(:,1) = PO4adsa_prev;
  PO4adsb = zeros(n,m);
  PO4adsb(:,1) = PO4adsb_prev;
  Ca2 = zeros(n,m);
  Ca2(:,1) = Ca2_prev;
  Ca3PO42 = zeros(n,m);
  Ca3PO42(:,1) = Ca3PO42_prev;
  OMS = zeros(n,m);
  OMS(:,1) = OMS_prev;
  H = zeros(n,m);
  H(:,1) = H_prev;
  OH = zeros(n,m);
  OH(:,1) = OH_prev;
  CO2 = zeros(n,m);
  CO2(:,1) = CO2_prev;
  CO3 = zeros(n,m);
  CO3(:,1) = CO3_prev;
  HCO3 = zeros(n,m);
  HCO3(:,1) = HCO3_prev;
  NH3 = zeros(n,m);
  NH3(:,1) = NH3_prev;
  NH4 = zeros(n,m);
  NH4(:,1) = NH4_prev;
  HS = zeros(n,m);
  HS(:,1) = HS_prev;
  H2S = zeros(n,m);
  H2S(:,1) = H2S_prev;
  H2CO3 = zeros(n,m);
  H2CO3(:,1) = H2CO3_prev;
  DOM1 = zeros(n,m);
  DOM1(:,1) = DOM1_prev;
  DOM2 = zeros(n,m);
  DOM2(:,1) = DOM2_prev;



  % Stiffness:
  % =======================================================================================================
  % CFL = (D_O2+Db)*dt/dx^2;


  % Solving equations!!!
  % =========================================================================================================

  for i=2:m
    C0 = [Ox(:,i-1), OM(:,i-1), OMb(:,i-1), NO3(:,i-1), FeOH3(:,i-1), SO4(:,i-1), NH4(:,i-1), Fe2(:,i-1), FeOOH(:,i-1), H2S(:,i-1), HS(:,i-1), FeS(:,i-1), S0(:,i-1), PO4(:,i-1), S8(:,i-1), FeS2(:,i-1), AlOH3(:,i-1), PO4adsa(:,i-1), PO4adsb(:,i-1), Ca2(:,i-1), Ca3PO42(:,i-1), OMS(:,i-1), H(:,i-1), OH(:,i-1), CO2(:,i-1), CO3(:,i-1), HCO3(:,i-1), NH3(:,i-1), H2CO3(:,i-1), DOM1(:,i-1), DOM2(:,i-1)];

      if any(any(isnan(C0)))
          error('NaN')
      end

      ts_during_one_dt = 1;
      int_method = 0;
      C_new = sediments_chemical_reactions_module(C0,dt,ts_during_one_dt, int_method);

      Ox(:,i-1)      = C_new(:,1);
      OM(:,i-1)      = C_new(:,2);
      OMb(:,i-1)     = C_new(:,3);
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
      DOM1(:,i-1)   = C_new(:,30);
      DOM2(:,i-1)   = C_new(:,31);



    % =======================================================================================================
    % Updating matrices and Solving eq-s
    % =======================================================================================================

      Ox(:,i) = pde_solver_solute(Ox_AL, Ox_AR, Ox(:,i-1), sediment_bc.Ox_c);
      OM(:,i) = pde_solver_solid(OM_AL, OM_AR, OM(:,i-1), sediment_bc.OM1_fx, sediment_params.solid_flux_coef);
      OMb(:,i) = pde_solver_solid(OMb_AL, OMb_AR, OMb(:,i-1), sediment_bc.OM2_fx, sediment_params.solid_flux_coef);
      NO3(:,i) = pde_solver_solute(NO3_AL, NO3_AR, NO3(:,i-1), sediment_bc.NO3_c);
      FeOH3(:,i) = pde_solver_solid(FeOH3_AL, FeOH3_AR, FeOH3(:,i-1), sediment_bc.FeOH3_fx, sediment_params.solid_flux_coef);
      SO4(:,i) = pde_solver_solute(SO4_AL, SO4_AR, SO4(:,i-1), sediment_bc.SO4_c);
      NH4(:,i) = pde_solver_solute(NH4_AL, NH4_AR, NH4(:,i-1), sediment_bc.NH4_c);
      Fe2(:,i) = pde_solver_solute(Fe2_AL, Fe2_AR, Fe2(:,i-1), sediment_bc.Fe2_c);
      FeOOH(:,i) = pde_solver_solid(FeOOH_AL, FeOOH_AR, FeOOH(:,i-1), sediment_bc.FeOOH_fx, sediment_params.solid_flux_coef);
      H2S(:,i) = pde_solver_solute(H2S_AL, H2S_AR, H2S(:,i-1), sediment_bc.H2S_c);
      HS(:,i) = pde_solver_solute(HS_AL, HS_AR, HS(:,i-1), sediment_bc.HS_c);
      FeS(:,i) = pde_solver_solid(FeS_AL, FeS_AR, FeS(:,i-1), sediment_bc.FeS_fx, sediment_params.solid_flux_coef);
      S0(:,i) = pde_solver_solute(S0_AL, S0_AR, S0(:,i-1), sediment_bc.S0_c);
      PO4(:,i) = pde_solver_solute(PO4_AL, PO4_AR, PO4(:,i-1), sediment_bc.PO4_c);
      S8(:,i) = pde_solver_solid(S8_AL, S8_AR, S8(:,i-1), sediment_bc.S8_fx, sediment_params.solid_flux_coef);
      FeS2(:,i) = pde_solver_solid(FeS2_AL, FeS2_AR, FeS2(:,i-1), sediment_bc.FeS2_fx, sediment_params.solid_flux_coef);
      AlOH3(:,i) = pde_solver_solid(AlOH3_AL, AlOH3_AR, AlOH3(:,i-1), sediment_bc.AlOH3_fx, sediment_params.solid_flux_coef);
      PO4adsa(:,i) = pde_solver_solid(PO4adsa_AL, PO4adsa_AR, PO4adsa(:,i-1), sediment_bc.PO4adsa_fx, sediment_params.solid_flux_coef);
      PO4adsb(:,i) = pde_solver_solid(PO4adsb_AL, PO4adsb_AR, PO4adsb(:,i-1), sediment_bc.PO4adsb_fx, sediment_params.solid_flux_coef);
      Ca2(:,i) = pde_solver_solute(Ca2_AL, Ca2_AR, Ca2(:,i-1), sediment_bc.Ca2_c);
      Ca3PO42(:,i) = pde_solver_solid(Ca3PO42_AL, Ca3PO42_AR, Ca3PO42(:,i-1), sediment_bc.Ca3PO42_fx, sediment_params.solid_flux_coef);
      OMS(:,i) = pde_solver_solid(OMS_AL, OMS_AR, OMS(:,i-1), sediment_bc.OMS_fx, sediment_params.solid_flux_coef);
      H(:,i) = pde_solver_solute(H_AL, H_AR, H(:,i-1), sediment_bc.H_c);
      OH(:,i) = pde_solver_solute(OH_AL, OH_AR, OH(:,i-1), sediment_bc.OH_c);
      CO2(:,i) = pde_solver_solute(CO2_AL, CO2_AR, CO2(:,i-1), sediment_bc.CO2_c);
      CO3(:,i) = pde_solver_solute(CO3_AL, CO3_AR, CO3(:,i-1), sediment_bc.CO3_c);
      HCO3(:,i) = pde_solver_solute(HCO3_AL, HCO3_AR, HCO3(:,i-1), sediment_bc.HCO3_c);
      NH3(:,i) = pde_solver_solute(NH3_AL, NH3_AR, NH3(:,i-1), sediment_bc.NH3_c);
      H2CO3(:,i) = pde_solver_solute(H2CO3_AL, H2CO3_AR, H2CO3(:,i-1), sediment_bc.H2CO3_c);
      DOM1(:,i) = pde_solver_solute(DOM1_AL, DOM1_AR, DOM1(:,i-1), sediment_bc.DOM1_c);
      DOM2(:,i) = pde_solver_solute(DOM2_AL, DOM2_AR, DOM2(:,i-1), sediment_bc.DOM2_c);



    % Add new species before this line.
    % =======================================================================================================


    % pH Module
    if sediment_params.pH_algorithm ~= 0
      [H(:,i), OH(:,i), DOM2(:,i), HCO3(:,i), CO2(:,i), CO3(:,i), NH3(:,i), NH4(:,i), HS(:,i), H2S(:,i)] = pH_module(sediment_params.pH_algorithm, H(:,i), OH(:,i), H2CO3(:,i), HCO3(:,i), CO2(:,i), CO3(:,i), NH3(:,i), NH4(:,i), HS(:,i), H2S(:,i), Fe2(:,i), Ca2(:,i), NO3(:,i), SO4(:,i), PO4(:,i), FeS(:,i), FeS2(:,i), FeOH3(:,i), FeOOH(:,i), Ca3PO42(:,i), PO4adsa(:,i), PO4adsb(:,i), sediment_bc.T);
    end

  end

% Estimate flux
  sediment_transport_fluxes.Ox           = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(Ox(:, end), D_O2, dx, fi), 31998);
  sediment_transport_fluxes.OM1          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( -sediment_bc.OM1_fx, 30973.762);
  sediment_transport_fluxes.OM2          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( -sediment_bc.OM2_fx, 12010.7);
  sediment_transport_fluxes.PO4          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(PO4(:, end), D_PO4, dx, fi), 30973.762);
  sediment_transport_fluxes.NO3          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(NO3(:, end), D_NO3, dx, fi), 62004);
  sediment_transport_fluxes.FeOH3        = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( -sediment_bc.FeOH3_fx, 106867.0);
  sediment_transport_fluxes.Fe2          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(Fe2(:, end), D_Fe2, dx, fi), 55845);
  sediment_transport_fluxes.NH4          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(NH4(:, end), D_NH4, dx, fi), 18038);
  sediment_transport_fluxes.AlOH3        = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( -sediment_bc.AlOH3_fx, 78003.6);
  sediment_transport_fluxes.PO4adsa      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( -sediment_bc.PO4adsa_fx, 30973.762);
  sediment_transport_fluxes.PO4adsb      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( -sediment_bc.PO4adsb_fx, 30973.762);
  sediment_transport_fluxes.SO4          = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(SO4(:, end), D_SO4, dx, fi), 96062);
  sediment_transport_fluxes.DOM1         = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(DOM1(:, end), D_DOM1, dx, fi), 30973.762);
  sediment_transport_fluxes.DOM2         = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( top_sediment_diffusion_flux(DOM2(:, end), D_DOM2, dx, fi), 12010.7);


  sediment_bioirrigation_fluxes.Ox   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(Ox(:, end), alfax, fi), x), 31998);
  sediment_bioirrigation_fluxes.PO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(PO4(:, end),  alfax,  fi), x), 30973.762);
  sediment_bioirrigation_fluxes.Fe2  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(Fe2(:, end),  alfax,  fi), x), 55845);
  sediment_bioirrigation_fluxes.NO3  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(NO3(:, end),  alfax,  fi), x), 62004);
  sediment_bioirrigation_fluxes.NH4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(NH4(:, end),  alfax,  fi), x), 18038);
  sediment_bioirrigation_fluxes.SO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(SO4(:, end),  alfax,  fi), x), 96062);
  sediment_bioirrigation_fluxes.DOM1 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(DOM1(:, end),  alfax,  fi), x), 30973.762);
  sediment_bioirrigation_fluxes.DOM2 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d( integrate_over_depth_2( bioirrigation(DOM2(:, end),  alfax,  fi), x), 12010.7);



  sediment_concentrations.Oxygen = Ox(:,end);
  sediment_concentrations.OM1 = OM(:,end);
  sediment_concentrations.OM2 = OMb(:,end);
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
  sediment_concentrations.DOM1 = DOM1(:,end);
  sediment_concentrations.DOM2 = DOM2(:,end);

  R_values = {
    % daily_average(R1_graphs, m) ,              'R1';
    % integrate_over_depth(R1_graphs, dx, m),    'R1 integrated';
    % daily_average(R2_graphs, m) ,              'R2';
    % integrate_over_depth(R2_graphs, dx, m),    'R2 integrated';
    % daily_average(R3_graphs,m),                'R3';
    % integrate_over_depth(R3_graphs, dx, m),    'R3 integrated';
    % daily_average(R4_graphs,m),                'R4';
    % integrate_over_depth(R4_graphs, dx, m),    'R4 integrated';
    % daily_average(R5_graphs, m),               'R5';
    % integrate_over_depth(R5_graphs, dx, m),    'R5 integrated';
    % daily_average(R6_graphs,m),                'R6';
    % integrate_over_depth(R6_graphs, dx, m),    'R6 integrated';
    };

  % Estimating of the water-column and sediment interaction due to diffusion and bioirrigation:


  sediment_additional_results = {...
    0,   'Ox integrated over depth flux to sediments';
  };




    if any(isnan(sediment_transport_fluxes.Ox))| any(isnan(sediment_bc.OM1_fx))| any(isnan(sediment_bc.OM2_fx))| any(isnan(sediment_bc.FeOH3_fx))| any(isnan(Ox)) | any(isnan(OM)) | any(isnan(OMb)) | any(isnan(NO3)) | any(isnan(FeOH3)) | any(isnan(SO4)) | any(isnan(NH4)) | any(isnan(Fe2)) | any(isnan(FeOOH)) | any(isnan(H2S)) | any(isnan(HS)) | any(isnan(FeS)) | any(isnan(S0)) | any(isnan(PO4)) | any(isnan(S8)) | any(isnan(FeS2)) | any(isnan(AlOH3)) | any(isnan(PO4adsa)) | any(isnan(PO4adsb)) | any(isnan(H)) | any(isnan(Ca2)) | any(isnan(Ca3PO42)) | any(isnan(OMS)) | any(isnan(OH)) | any(isnan(HCO3)) | any(isnan(CO2)) | any(isnan(CO3)) | any(isnan(NH3)) | any(isnan(H2CO3))
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

%% rk4: Runge-Kutta 4th order integration
function [C_new] = rk4(C0,ts, dt)
    % ts - how many time steps during 1 day
    dt = dt/ts;
    for i = 1:ts
        k_1 = dt.*sediment_rates(C0, dt);
        k_2 = dt.*sediment_rates(C0+0.5.*k_1, dt);
        k_3 = dt.*sediment_rates(C0+0.5.*k_2, dt);
        k_4 = dt.*sediment_rates(C0+k_3, dt);
        C_new = C0 + (k_1+2.*k_2+2.*k_3+k_4)/6;
        C0 = C_new;

        % if any(any(isnan(C_new)))
        %     error('NaN')
        % end
    end
end

%% butcher5: Butcher's Fifth-Order Runge-Kutta
function [C_new] = butcher5(C0,ts,dt)

    dt = dt/ts;
    for i = 1:ts
        k_1 = dt.*sediment_rates(C0, dt);
        k_2 = dt.*sediment_rates(C0 + 1/4.*k_1, dt);
        k_3 = dt.*sediment_rates(C0 + 1/8.*k_1 + 1/8.*k_2, dt);
        k_4 = dt.*sediment_rates(C0 - 1/2.*k_2 + k_3, dt);
        k_5 = dt.*sediment_rates(C0 + 3/16.*k_1 + 9/16.*k_4, dt);
        k_6 = dt.*sediment_rates(C0 - 3/7.*k_1 + 2/7.*k_2 + 12/7.*k_3 - 12/7.*k_4 + 8/7.*k_5, dt);
        C_new = C0 + (7.*k_1 + 32.*k_3 + 12.*k_4 + 32.*k_5 + 7.*k_6)/90;
        C0 = C_new;

        % if any(any(isnan(C_new)))
        %     error('NaN')
        % end
    end
end



function C_new = pde_solver_solute(AL, AR, C_old, const_bc)
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

function [C_new] = sediments_chemical_reactions_module(C0,dt,ts, method)
    % ts - how many time steps during 1 day
    if method == 0
        C_new = rk4(C0,ts,dt);
    elseif method == 1
        C_new = butcher5(C0,ts,dt);
    end
    C_new = (C_new>0).*C_new;
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

%% daily_average: returns the average rate during 1 day.
function [averaged] = daily_average(R)
  % R - rate of interest
  averaged = sum(R,2)/size(R,2);
end

function [flux] = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(fx, M_C)
  flux = fx * M_C * 10^4 / 365 / 10^6; % [umol/cm^2/y] -> [mg/m2/d]
end

function [flux] = top_sediment_diffusion_flux(C, D, dx, fi)
  % calculates flux of the particular dissolved specie through the top boundary of the sediment
  % in [ mg m-2 d-1 ] units
  % C(1) - BC of dissolved species
  % C - concentration
  % D - diffusion coefficient
  % M_C - molar mass in [ mg mol-1]
  % fi - porosity (no porosity because C is the concentration in pores (not bulk))

  % fourth-order
  flux = D * (-25 * C(1) + 48 * C(2) - 36 * C(3) + 16 * C(4) - 3 * C(5)) / dx / 12;  %  [umol/cm^2/y]

  % third order
  % flux = D * (-11 * C(1) + 18 * C(2) - 9 * C(3) + 2 * C(4)) / dx / 6;  %  [umol/cm^2/y]
  % flux = 0;  %  [umol/cm^2/y]

  % second order
  % flux = D * (-3 * C(1) + 4 * C(2) - C(3)) / dx / 2;  %  [umol/cm^2/y]

  % first order
  % flux = D * (C(3) - C(1)) / 2 / dx;  %  [umol/cm^2/y]
end

function bioR = bioirrigation(C, alfax, fi)
  % bioirrigation rate is the "artificial" function represents the bioirrigation by organisms in the sediment (worms etc) implemented according to Boudreau, B.P., 1999.
  % Co - bc wc-sediment value of current species
  % C - concentration profile of current species
  % fi - porosity
  Co = C(1);
  bioR = fi .* alfax .* (C - Co);
  % NOTE:Disabled?
  bioR = 0;
end


function [dcdt] = sediment_rates(C, dt)
% parameters for water-column chemistry
% NOTE: the rates are the same as in sediments. Units are per "year" due to time step is in year units too;

    global k_OM k_OMb k_DOM1 k_DOM2 Km_O2 Km_NO3 Km_FeOH3 Km_FeOOH Km_SO4 Km_oxao Km_amao Kin_O2 Kin_NO3  Kin_FeOH3 Kin_FeOOH k_amox k_Feox k_Sdis k_Spre k_FeS2pre k_pdesorb_c k_pdesorb_a k_pdesorb_b k_alum k_rhom   k_tS_Fe Ks_FeS k_Fe_dis k_Fe_pre k_apa  kapa k_oms k_tsox k_FeSpre f_pfe accel Cx1 Ny1 Pz1 Cx2 Ny2 Pz2 F Ny1 Ny2 Pz1 Pz2 alfax fi

    dcdt=zeros(size(C));

    if any(isnan(C))
      error('Breaking out of Sediments function: NaN values');
    end

    Ox      = C(:,1) .* (C(:,1)>0) ;
    OM      = C(:,2) .* (C(:,2)>0) ;
    OMb     = C(:,3) .* (C(:,3)>0) ;
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
    DOM1    = C(:,30) .* (C(:,30)>0) ;
    DOM2    = C(:,31) .* (C(:,31)>0) ;

    f_O2    = Ox ./  (Km_O2 + Ox);
    f_NO3   = NO3 ./  (Km_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    f_FeOH3 = FeOH3 ./  (Km_FeOH3 + FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    f_FeOOH = FeOOH ./  (Km_FeOOH + FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    f_SO4   = SO4 ./ (Km_SO4 + SO4 ) .* Kin_FeOOH ./ (Kin_FeOOH + FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + Ox);
    Sum_H2S = H2S + HS;;
    Sat_FeS = Fe2*1e-3 .* Sum_H2S*1e-3 ./ ((H*1e-3).^2 .* Ks_FeS);


    R1a = k_OM.*OM .* f_O2 * accel;
    R1b = k_OMb.*OMb .* f_O2 * accel;
    R1c = k_DOM1.* DOM1 .* f_O2 * accel;
    R1d = k_DOM2.*DOM2 .* f_O2 * accel;

    R2a = k_OM.*OM .* f_NO3 * accel;
    R2b = k_OMb.*OMb .* f_NO3 * accel;
    R2c = k_DOM1.*DOM1 .* f_NO3 * accel;
    R2d = k_DOM2.*DOM2 .* f_NO3 * accel;

    R3a = k_OM.*OM .* f_FeOH3;
    R3b = k_OMb .*OMb .* f_FeOH3;
    R3c = k_DOM1 .*DOM1 .* f_FeOH3;
    R3d = k_DOM2 .*DOM2 .* f_FeOH3;

    R4a = k_OM.*OM .* f_FeOOH;
    R4b = k_OMb.*OMb .* f_FeOOH;
    R4c = k_DOM1.*DOM1 .* f_FeOOH;
    R4d = k_DOM2.*DOM2 .* f_FeOOH;

    R5a = k_OM.*OM .* f_SO4  ;
    R5b = k_OMb.*OMb .* f_SO4 ;
    R5c = k_DOM1.*DOM1 .* f_SO4 ;
    R5d = k_DOM2.*DOM2 .* f_SO4 ;

    Ra = R1a+R2a+R3a+R4a+R5a;
    Rb = R1b+R2b+R3b+R4b+R5b;
    Rc = R1c+R2c+R3c+R4c+R5c;
    Rd = R1d+R2d+R3d+R4d+R5d;
    R1 = R1a+R1b+R1c+R1d;
    R2 = R2a+R2b+R2c+R2d;
    R3 = R3a+R3b+R3c+R3d;
    R4 = R4a+R4b+R4c+R4d;
    R5 = R5a+R5b+R5c+R5d;

    R6 = k_tsox * Ox .* Sum_H2S;
    R7 = k_tS_Fe * FeOH3 .*  Sum_H2S;

    R8 = k_Feox .* Fe2 .* Ox;
    % NOTE: Due to the reaction is too fast and could cause overshooting:
    % we need to make this check if R*dt > Conc of source:
    % if R*dt > Conc then R8 = C/dt
    % if R*dt < Conc then R8 = R8
    % R8 = (R8.*dt < Fe2/50).*R8 + (R8.*dt > Fe2/50).* R8 ./ 1000;
    R8 = (R8.*dt < Fe2).*R8 + (R8.*dt > Fe2).* Fe2 ./ (dt) * 0.5;
    R8 = (R8.*dt < Ox).*R8 + (R8.*dt > Ox).* Ox ./ (dt) * 0.5;

    % R9 = k_amox * Ox ./ (Km_oxao + Ox) .* (NH4 ./ (Km_amao + NH4)); % NOTE: Doesnt work - Highly unstable.
    R9 = k_amox  .* NH4 .* Ox;
    R9 = (R9.*dt < NH4).*R9 + (R9.*dt > NH4).* NH4 ./ (dt) * 0.5;
    R9 = (R9.*dt < Ox).*R9 + (R9.*dt > Ox).* Ox ./ (dt) * 0.5;

    R10a = k_oms * Sum_H2S .* OM;
    R10b = k_oms * Sum_H2S .* OMb;
    R10c = k_oms * Sum_H2S .* DOM1;
    R10d = k_oms * Sum_H2S .* DOM2;
    R11 = k_FeSpre .* FeS .* S0;
    R12 = k_rhom * Ox .* FeS;
    R13 = k_FeS2pre .* FeS .* Sum_H2S;
    R14a = k_Fe_pre .* (Sat_FeS - 1);
    R14b  =  k_Fe_dis .* FeS .* (1 - Sat_FeS);
    R14a = (R14a >= 0) .* R14a; % can only be non negative
    R14b  = (R14b >= 0) .* R14b; % can only be non negative
    R15a = k_Spre * S0;
    R15b = k_Sdis .* S8;

    R16a = k_pdesorb_a * (FeOH3 - PO4adsa) .* PO4;
    R16b = f_pfe .* (4 * R3 + 2 * R7);
    R16b = (R16b.*dt < PO4adsa).*R16b + (R16b.*dt > PO4adsa).* PO4adsa ./ (dt) * 0.5;
    R17a = k_pdesorb_b * (FeOOH - PO4adsb) .* PO4;
    R17b = f_pfe .* (4 * R4);
    R17b = (R17b.*dt < PO4adsb).*R17b + (R17b.*dt > PO4adsb).* PO4adsb ./ (dt) * 0.5;

    R18a = k_pdesorb_c .* PO4 .* AlOH3;
    R18b = 0;
    R19 = k_apa * (PO4 - kapa);
    R19 = (R19 >= 0) .* R19; % can only be non negative

    % for stoichiometry check:
    % Canavan, R. W., Slomp, C. P., Jourabchi, P., Van Cappellen, P., Laverman, A. M., & van den Berg, G. A. (2006). Organic matter mineralization in sediment of a coastal freshwater lake and response to salinization. Geochimica Et Cosmochimica Acta, 70(11), 2836–2855. http://doi.org/10.1016/j.gca.2006.03.012

    % F = 2*(1-fi)./fi;

    dcdt(:,1)  = -0.25 * R8  - 2 * R9  - (Cx1*R1a + Cx2*R1b) .* F - (Cx1*R1c + Cx2*R1d) - 3 * R12 + bioirrigation(Ox, alfax, fi); % Ox
    dcdt(:,2)  = -1*Ra - R10a; % POC1
    dcdt(:,3)  = -1*Rb - R10b; % POC2
    dcdt(:,4)  = - 0.8*(Cx1*R2a+Cx1*R2b) .* F - 0.8*(Cx1*R2c+Cx1*R2d)+ R9 + bioirrigation(NO3, alfax, fi); % NO3
    dcdt(:,5)  = -4 * (Cx1*R3a + Cx2*R3b + Cx1*R3c + Cx2*R3d) - 2*R7 + R8./ F - R16a; % FeOH3
    dcdt(:,6)  = - 0.5*(Cx1*R5a + Cx2*R5b) .* F -0.5*(Cx1*R5c + Cx2*R5d)+ R6 + bioirrigation(SO4, alfax, fi); % SO4
    dcdt(:,7)  = (Ny1 * Ra + Ny2 * Rb) .* F + (Ny1 * Rc + Ny2 * Rd) - R9 + bioirrigation(NH4, alfax, fi); % NH4
    dcdt(:,8)  = 4*(Cx1*R3a + Cx2*R3b) .* F + 4* (Cx1*R3c + Cx2*R3d) + 4*(Cx1*R4a + Cx2*R4b) .* F + 4 * (Cx1*R4c + Cx2*R4d) + 2*R7 - R8 + R14b - R14a + bioirrigation(Fe2, alfax, fi); % Fe2
    dcdt(:,9)  = -4*(Cx1*R4a + Cx2*R4b + Cx1*R4c + Cx2*R4d) + R12 - R17a; % FeOOH
    dcdt(:,10) = + bioirrigation(H2S, alfax, fi); % H2S
    dcdt(:,11) = 0.5*(Cx1*R5a + Cx2*R5b) .* F + 0.5 * (Cx1*R5c + Cx2*R5d) - R6 - R7 + R14b - R14a - R10a - R10b - R10c - R10d -R13 + bioirrigation(HS, alfax, fi); % HS
    dcdt(:,12) = - R14b - R11 - 4*R12 -R13 + R14a; % FeS
    dcdt(:,13) = - R11 - R15a + R7 + R15b; % S0
    dcdt(:,14) = (Pz1 * Ra + Pz2 * Rb) .* F + (Pz1 * Rc + Pz2 * Rd) + R16b + R17b - 2 * R19 - R18a - R16a - R17a + bioirrigation(PO4, alfax, fi); % PO4
    dcdt(:,15) = 4*R12 - R15b + R15a; % S8
    dcdt(:,16) = + R11 + R13; % FeS2
    dcdt(:,17) = -R18a; % AlOH3
    dcdt(:,18) = R16a - R16b; % PO4adsa
    dcdt(:,19) = R17a - R17b; % PO4adsb
    dcdt(:,20) = -3*R19 + bioirrigation(Ca2, alfax, fi); % Ca2
    dcdt(:,21) = R19; % Ca3PO42
    dcdt(:,22) = R10a + R10b + R10c + R10d; % OMS
    dcdt(:,23) = 0; % H
    dcdt(:,24) = 0; % OH
    dcdt(:,25) = ((Cx1 - Ny1 + 2*Pz1)*R1a + (Cx2 - Ny2 + 2*Pz2)*R1b  + (0.2*Cx1 - Ny1 + 2*Pz1)*R2a +  (0.2*Cx2 - Ny2 + 2*Pz2)*R2b - (7*Cx1 + Ny1 + 2*Pz1)*(R3a+R4a) - (7*Cx2 + Ny2 + 2*Pz2)*(R3b+R4b)  - (Ny1 - 2*Pz1)*R5a + (Ny2 - 2*Pz2)*R5b) .* F  +  (Cx1 - Ny1 + 2*Pz1)*R1c + (Cx2 - Ny2 + 2*Pz2)*R1d + (0.2*Cx1 - Ny1 + 2*Pz1)*R2c +  (0.2*Cx2 - Ny2 + 2*Pz2)*R2d - (7*Cx1 + Ny1 + 2*Pz1)*(R3c+R4c) - (7*Cx2 + Ny2 + 2*Pz2)*(R3d+R4d)  - (Ny1 - 2*Pz1)*R5c + (Ny2 - 2*Pz2)*R5d + 2*R8 + 2*R9 + bioirrigation(CO2, alfax, fi);  % CO2
    dcdt(:,26) + bioirrigation(CO3, alfax, fi); % CO3
    dcdt(:,27) = ((0.8*Cx1 + Ny1 - 2*Pz1)*R2a + (0.8*Cx2 + Ny2 - 2*Pz2)*R2b  + (8*Cx1+Ny1-2*Pz1)*(R3a + R4a) +(8*Cx2+Ny2-2*Pz2)*(R3b + R4b)  + (Cx1+Ny1-2*Pz1)*R5a + (1*Cx2+Ny2-2*Pz2)*R5b ) .* F + (0.8*Cx1 + Ny1 - 2*Pz1)*R2c + (0.8*Cx2 + Ny2 - 2*Pz2)*R2d + (8*Cx1+Ny1-2*Pz1)*(R3c + R4c) +(8*Cx2+Ny2-2*Pz2)*(R3d + R4d) + (Cx1+Ny1-2*Pz1)*R5c + (1*Cx2+Ny2-2*Pz2)*R5d -  2*R8 - 2*R9 + bioirrigation(HCO3, alfax, fi); % HCO3

    dcdt(:,28) + bioirrigation(NH3, alfax, fi); % NH3
    dcdt(:,29) + bioirrigation(H2CO3, alfax, fi); % H2CO3
    dcdt(:,30) = -1*Rc - R10c + bioirrigation(DOM1, alfax, fi); % DOM1
    dcdt(:,31) = -1*Rd - R10d + bioirrigation(DOM2, alfax, fi); % DOM2
end

