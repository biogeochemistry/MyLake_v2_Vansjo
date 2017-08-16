% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2005
%
% VERSION 1.2.1 (two phytoplankton groups are included; variable Cz denotes
% this second group now. Frazil ice included + some small bug-fixes and code rearrangements. Using convection_v2.m code)
%
% Main module
% Code checked by TSA, xx.xx.200x
% Last modified by TSA, 21.08.2007

% Modified to include Fokema-module by Kai Rasmus. 16.5.2007
% Modified to include the latest Fokema module 30.12.2010 by PK
% New matrices: DOCzt1,DOCzt2,DOCzt3,Daily_BB1t,Daily_BB2t,Daily_BB3t,Daily_PBt

% New DIC variable 29.12.2010 (incl. inflow, convection, diffusion) by PK
% New O2 variable 10.2.2011 by PK

function [MyLake_results, sediment_results] = ...
    solvemodel_v2(M_start,M_stop,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,varargin)

warning off MATLAB:fzero:UndeterminedSyntax %suppressing a warning message
% Inputs (to function)
%       M_start     : Model start date [year, month, day]
%       M_stop      : Model stop date [year, month, day]
%               + Input filenames and sheetnames

% Inputs (received from input module):
%		tt		: Solution time domain (day)
%       In_Z    : Depths read from initial profiles file (m)
%       In_Az   : Areas read from initial profiles file (m2)
%       In_Tz   : Initial temperature profile read from initial profiles file (deg C)
%       In_Cz   : Initial chlorophyll (group 2) profile read from initial profiles file (-)
%       In_POCz   : Initial sedimenting tracer (or suspended inorganic matter) profile read from initial profiles file (kg m-3)
%       In_TPz  : Initial total P profile read from initial profiles file (incl. DOP & Chla & Cz) (mg m-3)
%       In_DOPz  : Initial dissolved organic P profile read from initial profiles file (mg m-3)
%       In_Chlz : Initial chlorophyll (group 1) profile read from initial profiles file (mg m-3)
%       In_DOCz  : Initial DOC profile read from initial profiles file (mg m-3)
%       In_DICz  : Initial DIC profile read from initial profiles file (mg m-3) (PK)
%       In_O2z  : Initial oxygen profile read from initial profiles file (mg m-3) (PK)
%       In_TPz_sed  : Initial total P profile in the sediment compartments read from initial profiles file (mg m-3)
%       In_Chlz_sed : Initial chlorophyll profile  (groups 1+2) in the sediment compartments read from initial profiles file (mg m-3)
%       In_FIM      : Initial profile of volume fraction of inorganic matter in the sediment solids (dry weight basis)
%       Ice0            : Initial conditions, ice and snow thicknesses (m) (Ice, Snow)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 23 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (23 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 23 parameters that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (23 * 2)
%       Bio_par_names   : Names for Bio_par

% Outputs (other than Inputs from input module):
%		Qst : Estimated surface heat fluxes ([sw, lw, sl] * tt) (W m-2)
%		Kzt	: Predicted vertical diffusion coefficient (tt * zz) (m2 d-1)
%		Tzt	: Predicted temperature profile (tt * zz) (deg C)
%		Czt	: Predicted chlorophyll (group 2) profile (tt * zz) (-)
%		POCzt	: Predicted passive sedimenting tracer (or suspended inorganic matter) profile (tt * zz) (kg m-3)=(g L-1)
%		Pzt	: Predicted dissolved inorganic phosphorus profile (tt * zz) (mg m-3)
%		Chlzt	    : Predicted chlorophyll (group 1) profileo (tt * zz) (mg m-3)
%		PPzt	    : Predicted particulate inorganic phosphorus profile (tt * zz) (mg m-3)
%		DOPzt	    : Predicted dissolved organic phosphorus profile (tt * zz) (mg m-3)
%		DOCzt	    : Predicted dissolved organic carbon (DOC) profile (tt * zz) (mg m-3)
%		DICzt	    : Predicted dissolved inorganic carbon (DIC) profile (tt * zz) (mg m-3) (PK)
%		CO2zt	    : Predicted dissolved carbon dioxide profile (tt * zz) (mg m-3) (PK)
%		O2zt	    : Predicted dissolved oxygen profile (tt * zz) (mg m-3) (PK)
%       O2_sat_rel  : Predicted relative oxygen saturation (PK)
%       O2_sat_abs  : Predicted absolute oxygen saturation (PK)
%		Qz_sed      : Predicted  sediment-water heat flux (tt * zz) (W m-2, normalised to lake surface area)
%       lambdazt    : Predicted average total light attenuation coefficient down to depth z (tt * zz) (m-1)
%       P3zt_sed    : Predicted P conc. in sediment for P (mg m-3), PP(mg kg-1 dry w.) and Chl (mg kg-1 dry w.) (tt * zz * 3)
%       P3zt_sed_sc : Predicted P source from sediment for P, PP and Chl (mg m-3 day-1) (tt * zz * 3)
%       His         : Ice information matrix ([Hi Hs Hsi Tice Tair rho_snow IceIndicator] * tt)
%       DoF, DoM    : Days of freezing and melting (model timestep number)
%       MixStat     : Temporary variables used in model testing, see code (N * tt)

% Fokema outputs
%       CDOMzt      : Coloured dissolved organic matter absorption m-1
%                   : (tt * zz)

% These variables are still global and not transferred by functions
global ies80


tic
disp(['Running MyLake-Sediment from ' datestr(datenum(M_start)) ' to ' datestr(datenum(M_stop)) ' ...']);

% ===Switches===
snow_compaction_switch=1;       %snow compaction: 0=no, 1=yes
river_inflow_switch=1;          %river inflow: 0=no, 1=yes
deposition_switch= 0;			%human impact, atm deposition , point source addition %% NEW_DOCOMO
sediment_heatflux_switch=1;     %heatflux from sediments: 0=no, 1=yes
selfshading_switch=1;           %light attenuation by chlorophyll a: 0=no, 1=yes
tracer_switch=1;                %simulate tracers:  0=no, 1=yes
matsedlab_sediment_module = 1;  %MATSEDLAB sediment module  %% NEW_DOCOMO
wc_chemistry_module = 1;        % WC chemistry module: on/off
wc_int_method = 0;              % WC chemistry module: method: 0 = Runge-Kutta 4th order; 1 = Buthcer's 5th Order;
%fokema
photobleaching=0;               %photo bleaching: 0=TSA model, 1=FOKEMA model
floculation_switch=1;           % floculation according to Wachenfeldt 2008  %% NEW_DOCOMO
resuspension_enabled=0;         % Resuspension switch

rate_estimator_switch=0;        % estimate rates or not, additional cost of about 20% of computational time;
% ==============

dt=1.0; %model time step = 1 day (DO NOT CHANGE!)

if (nargin>8) %if optional command line parameter input is used
    disp('Bypassing input files...Running with input data & parameters given on command line');
    [In_Z,In_Az,tt,In_Tz,In_Cz,In_POCz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_DICz,In_O2z,In_NO3z,In_NH4z,In_SO4z,In_HSz,In_H2Sz,In_Fe2z,In_Ca2z,In_pHz,In_CH4z,In_Fe3z,In_Al3z,In_SiO4z,In_SiO2z,In_diatomz,In_POPz,In_TPz_sed,In_Chlz_sed,In_FIM,Ice0,Wt,Inflw,...
        Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names, Deposition]...
        = deal(varargin{:});
else
    %Read input data
    [In_Z,In_Az,tt,In_Tz,In_Cz,In_POCz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_DICz,In_TPz_sed,In_Chlz_sed,In_O2z,In_NO3z,In_NH4z,In_SO4z,In_HSz,In_H2Sz,In_Fe2z,In_Ca2z,In_pHz,In_CH4z,In_Fe3z,In_Al3z,In_SiO4z,In_SiO2z,In_diatomz,In_POPz,In_FIM,Ice0,Wt,Inflw,...
         Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names]...
        = modelinputs_v2(M_start,M_stop,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,dt);
end

load albedot1.mat; %load albedot1 table, in order to save execution time

% Unpack the more fixed parameter values from input array "Phys_par"
dz = Phys_par(1); %grid step size (m)

zm = In_Z(end); %max depth
zz = [0:dz:zm-dz]'; %solution depth domain

Kz_K1 = Phys_par(2); % open water diffusion parameter (-)
Kz_K1_ice = Phys_par(3); % under ice diffusion parameter (-)
Kz_N0 = Phys_par(4); % min. stability frequency (s-2)
C_shelter = Phys_par(5); % wind shelter parameter (-)
lat = Phys_par(6); %latitude (decimal degrees)
lon = Phys_par(7); %longitude (decimal degrees)
alb_melt_ice = Phys_par(8);   %albedo of melting ice (-)
alb_melt_snow = Phys_par(9); %albedo of melting snow (-)
PAR_sat = Phys_par(10);         %PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
f_par = Phys_par(11);           %Fraction of PAR in incoming solar radiation (-)
beta_chl = Phys_par(12);        %Optical cross_section of chlorophyll (m2 mg-1)
lambda_i = Phys_par(13);       %PAR light attenuation coefficient for ice (m-1)
lambda_s = Phys_par(14);       %PAR light attenuation coefficient for snow (m-1)
F_sed_sld = Phys_par(15);      %volume fraction of solids in sediment (= 1-porosity)
I_scV = Phys_par(16); %scaling factor for inflow volume (-)
I_scT = Phys_par(17); %scaling coefficient for inflow temperature (-)
I_scC = Phys_par(18); %scaling factor for inflow concentration of C (-)
I_scPOC = Phys_par(19); %scaling factor for inflow concentration of POC (-)
I_scTP = Phys_par(20); %scaling factor for inflow concentration of total P (-)
I_scDOP = Phys_par(21); %scaling factor for inflow concentration of diss. organic P (-)
I_scChl = Phys_par(22); %scaling factor for inflow concentration of Chl a (-)
I_scDOC = Phys_par(23); %scaling factor for inflow concentration of DOC  (-)
I_scPOP = Phys_par(24); %scaling factor for inflow concentration of POP  (-)
I_scO = Phys_par(25); %scaling factor for inflow concentration of O2  (-)
I_scDIC = Phys_par(26);   %Scaling factor for inflow concentration of DIC  (-)
I_scNO3 = Phys_par(27);   %scaling factor for inflow concentration   (-)
I_scNH4 = Phys_par(28);    %scaling factor for inflow concentration   (-)
I_scSO4 = Phys_par(29);    %scaling factor for inflow concentration   (-)
I_scFe2 = Phys_par(30);    %scaling factor for inflow concentration   (-)
I_scCa2 = Phys_par(31);    %scaling factor for inflow concentration   (-)
I_scpH = Phys_par(32);    %scaling factor for inflow concentration   (-)
I_scCH4 = Phys_par(33);    %scaling factor for inflow concentration   (-)
I_scFe3 = Phys_par(34);    %scaling factor for inflow concentration   (-)
I_scAl3 = Phys_par(35);    %scaling factor for inflow concentration   (-)
I_scSiO4 = Phys_par(36);    %scaling factor for inflow concentration   (-)
I_scSiO2 = Phys_par(37);    %scaling factor for inflow concentration   (-)
I_scdiatom = Phys_par(38);    %scaling factor for inflow concentration   (-)

% Unpack the more site specific parameter values from input array "Bio_par"

swa_b0 = Bio_par(1); % non-PAR light attenuation coeff. (m-1)
swa_b1 = Bio_par(2); %  PAR light attenuation coeff. (m-1)
S_res_epi = Bio_par(3);      %Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
S_res_hypo = Bio_par(4);     %Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
H_sed = Bio_par(5);          %height of active sediment layer (m, wet mass)
Psat_L = Bio_par(6);           %Half saturation parameter for Langmuir isotherm
Fmax_L = Bio_par(7);    %Scaling parameter for Langmuir isotherm !!!!!!!!!!!!

w_s = Bio_par(8);              %settling velocity for S (m day-1)
w_chl = Bio_par(9);            %settling velocity for Chl a (m day-1)
Y_cp = Bio_par(10);            %yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)
m_twty = Bio_par(11);          %loss rate (1/day) at 20 deg C
g_twty = Bio_par(12);          %specific growth rate (1/day) at 20 deg C
k_twty = Bio_par(13);          %specific Chl a to P transformation rate (1/day) at 20 deg C
dop_twty = Bio_par(14);        %specific DOP to P transformation rate (day-1) at 20 deg C
P_half = Bio_par(15);          %Half saturation growth P level (mg/m3)

%NEW!!!===parameters for the 2 group of chlorophyll variable
PAR_sat_2 = Bio_par(16);        %PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
beta_chl_2 = Bio_par(17);       %Optical cross_section of chlorophyll (m2 mg-1)
w_chl_2 = Bio_par(18);          %Settling velocity for Chl a (m day-1)
m_twty_2 = Bio_par(19);         %Loss rate (1/day) at 20 deg C
g_twty_2 = Bio_par(20);         %Specific growth rate (1/day) at 20 deg C
P_half_2 = Bio_par(21);         %Half saturation growth P level (mg/m3)

oc_DOC = Bio_par(22);           %Optical cross-section of DOC (m2/mg DOC)
qy_DOC = Bio_par(23);           %Quantum yield (mg DOC degraded/mol quanta)
%===========

% Parameters for oxygen

k_BOD = Bio_par(24);            %Organic decomposition rate (1/d)
k_SOD = Bio_par(25);            %Sedimentary oxygen demand (mg m-2 d-1)
theta_bod = Bio_par(26);        %Temperature adjustment coefficient for BOD, T ? 10 °C
theta_bod_ice = Bio_par(27);    %Temperature adjustment coefficient for BOD, T < 10 °C
theta_sod = Bio_par(28);        %Temperature adjustment coefficient for SOD, T ? 10 °C
theta_sod_ice = Bio_par(29);    %Temperature adjustment coefficient for SOD, T < 10 °C
BOD_temp_switch = Bio_par(30);             %Threshold for bod or bod_ice °C
pH = Bio_par(31);               %Lake water pH
% WC chemistry:
Q10 = Bio_par(32);
wc_factor = Bio_par(33);
T_ref = Bio_par(34);



% ====== Other variables/parameters not read from the input file:

Nz=length(zz); %total number of layers in the water column
N_sed=26; %total number of layers in the sediment column

theta_m = exp(0.1*log(2));    %loss and growth rate parameter base, ~1.072
e_par = 240800;               %Average energy of PAR photons (J mol-1)

% diffusion parameterization exponents
Kz_b1 = 0.43;
Kz_b1_ice  = 0.43;

% ice & snow parameter values
rho_fw=1000;        %density of freshwater (kg m-3)
rho_ice=910;        %ice (incl. snow ice) density (kg m-3)
rho_new_snow=250;   %new-snow density (kg m-3)
max_rho_snow=450;   %maximum snow density (kg m-3)
L_ice=333500;       %latent heat of freezing (J kg-1)
K_ice=2.1;          %ice heat conduction coefficient (W m-1 K-1)
C1=7.0;             %snow compaction coefficient #1
C2=21.0;            %snow compaction coefficient #2

Tf=0;               %water freezing point temperature (deg C)

F_OM=1e+6*0.012;    %mass fraction [mg kg-1] of P of dry organic matter (assuming 50% of C, and Redfield ratio)

K_sed=0.035;      %thermal diffusivity of the sediment (m2 day-1)
rho_sed=2500;      %bulk density of the inorganic solids in sediment (kg m-3)
rho_org=1000;      %bulk density of the organic solids in sediment (kg m-3)
cp_sed=1000;       %specific heat capasity of the sediment (J kg-1 K-1)

ksw=1e-3; %sediment pore water mass transfer coefficient (m/d)
Fmax_L_sed=Fmax_L;
Fstable=655; % Inactive P conc. in inorg. particles (mg/kg dw);

Frazil2Ice_tresh=0.03;  % treshold (m) where frazil is assumed to turn into a solid ice cover NEW!!!


% Allocate and initialise output data matrices
Qst = zeros(3,length(tt));
Kzt = zeros(Nz,length(tt));
Tzt = zeros(Nz,length(tt));
Czt = zeros(Nz,length(tt));
POCzt = zeros(Nz,length(tt));
Pzt = zeros(Nz,length(tt));
Chlzt = zeros(Nz,length(tt));
PPzt = zeros(Nz,length(tt));
DOPzt = zeros(Nz,length(tt));
DOCzt = zeros(Nz,length(tt));
DICzt = zeros(Nz,length(tt));
CO2zt = zeros(Nz,length(tt));
O2zt = zeros(Nz,length(tt));

NO3zt  = zeros(Nz,length(tt));
NH4zt  = zeros(Nz,length(tt));
SO4zt  = zeros(Nz,length(tt));
HSzt  = zeros(Nz,length(tt));
H2Szt  = zeros(Nz,length(tt));
Fe2zt  = zeros(Nz,length(tt));
Ca2zt  = zeros(Nz,length(tt));
pHzt  = zeros(Nz,length(tt));
CH4zt  = zeros(Nz,length(tt));
Fe3zt  = zeros(Nz,length(tt));
Al3zt  = zeros(Nz,length(tt));
SiO4zt  = zeros(Nz,length(tt));
SiO2zt  = zeros(Nz,length(tt));
diatomzt  = zeros(Nz,length(tt));
POPzt  = zeros(Nz,length(tt));


H_sw_zt  = zeros(Nz,length(tt));
H_sw_zt_2  = zeros(Nz,length(tt));
PAR_zt  = zeros(Nz,length(tt));

O2_diffzt = zeros(Nz,length(tt));
O2_sat_relt = zeros(Nz,length(tt));
O2_sat_abst = zeros(Nz,length(tt));
Qzt_sed = zeros(Nz,length(tt));
lambdazt = zeros(Nz,length(tt));
P3zt_sed = zeros(Nz,length(tt),4); %3-D
P3zt_sed_sc = zeros(Nz,length(tt),3); %3-D
His = zeros(8,length(tt)); %NEW!!!
MixStat = zeros(23,length(tt));
% Fokema
CDOMzt=zeros(Nz,length(tt));
DOCzt1=zeros(Nz,length(tt)); %Fokema-model subpool 1
DOCzt2=zeros(Nz,length(tt)); %Fokema-model subpool 2
DOCzt3=zeros(Nz,length(tt)); %Fokema-model subpool 3
DOC1tfrac=zeros(Nz,length(tt)); %Fokema-model subpool 1
DOC2tfrac=zeros(Nz,length(tt)); %Fokema-model subpool 2 fraction
DOC3tfrac=zeros(Nz,length(tt)); %Fokema-model subpool 3 fraction
Daily_BB1t=zeros(Nz,length(tt)); %Fokema-model subpool 1 daily bacterial decomposition
Daily_BB2t=zeros(Nz,length(tt)); %Fokema-model subpool 2 daily bacterial decomposition
Daily_BB3t=zeros(Nz,length(tt)); %Fokema-model subpool 3 daily bacterial decomposition
Daily_PBt=zeros(Nz,length(tt)); %Fokema-model daily photobleaching

surfaceflux = zeros(1,length(tt)); %CO2 surface flux
CO2_eqt = zeros(1,length(tt));     %CO2 equilibrium concentration
CO2_ppmt = zeros(1,length(tt));    %CO2 fraction in air
K0t = zeros(1,length(tt));         %CO2 solubility coefficient

O2fluxt = zeros(1,length(tt));     %oxygen surface flux
O2_eqt = zeros(1,length(tt));      %O2 equilibrium concentration
K0_O2t = zeros(1,length(tt));      %O2 solubility coefficient
dO2Chlt = zeros(Nz,length(tt));    %Oxygen change due to phytoplankton (mg m-3))
dO2BODt = zeros(Nz,length(tt));    %Oxygen consumption due to BOD (mg m-3))
% dO2SODt = zeros(Nz,length(tt));    %Oxygen consumption due to SOD (mg m-3))
dfloc_DOC =  zeros(Nz,length(tt));  % floculation rates
testi1t = zeros(Nz,length(tt));
testi2t = zeros(Nz,length(tt));testi3t = zeros(Nz,length(tt));
lvlDzt = zeros(1,length(tt));

% Initial profiles

Az = interp1(In_Z,In_Az,zz);
Vz = dz * (Az + [Az(2:end); 0]) / 2;

T0 = interp1(In_Z,In_Tz,zz+dz/2); % Initial temperature distribution (deg C)
C0 = interp1(In_Z,In_Cz,zz+dz/2); % Initial  chlorophyll (group 2) distribution (mg m-3)
POC0 = interp1(In_Z,In_POCz,zz+dz/2); % Initial passive sedimenting tracer (or suspended inorganic matter) distribution (kg m-3)
TP0 = interp1(In_Z,In_TPz,zz+dz/2);	% Initial total P distribution (incl. DOP & Chla & Cz) (mg m-3)
DOP0 = interp1(In_Z,In_DOPz,zz+dz/2);	% Initial dissolved organic P distribution (mg m-3)
Chl0 = interp1(In_Z,In_Chlz,zz+dz/2);	% Initial chlorophyll (group 2) distribution (mg m-3)
DOC0 = interp1(In_Z,In_DOCz,zz+dz/2);	% Initial DOC distribution (mg m-3)
DIC0 = interp1(In_Z,In_DICz,zz+dz/2);   % Initial DIC distribution (mg m-3)
O20 = interp1(In_Z,In_O2z,zz+dz/2);   % Initial oxygen distribution (mg m-3)
TP0_sed = interp1(In_Z,In_TPz_sed,zz+dz/2); % Initial total P distribution in bulk wet sediment ((mg m-3); particles + porewater)
Chl0_sed = interp1(In_Z,In_Chlz_sed,zz+dz/2); % Initial chlorophyll (group 1+2) distribution in bulk wet sediment (mg m-3)
FIM0 = interp1(In_Z,In_FIM,zz+dz/2);     % Initial sediment solids volume fraction of inorganic matter (-)

NO30 = interp1(In_Z,In_NO3z,zz+dz/2);
NH40 = interp1(In_Z,In_NH4z,zz+dz/2);
SO40 = interp1(In_Z,In_SO4z,zz+dz/2);
HS0 = interp1(In_Z,In_HSz,zz+dz/2);
H2S0 = interp1(In_Z,In_H2Sz,zz+dz/2);
Fe20 = interp1(In_Z,In_Fe2z,zz+dz/2);
Ca20 = interp1(In_Z,In_Ca2z,zz+dz/2);
pH0 = interp1(In_Z,In_pHz,zz+dz/2);
CH40 = interp1(In_Z,In_CH4z,zz+dz/2);
Fe30 = interp1(In_Z,In_Fe3z,zz+dz/2);
Al30 = interp1(In_Z,In_Al3z,zz+dz/2);
SiO40 = interp1(In_Z,In_SiO4z,zz+dz/2);
SiO20 = interp1(In_Z,In_SiO2z,zz+dz/2);
diatom0 = interp1(In_Z,In_diatomz,zz+dz/2);
POP0 = interp1(In_Z,In_POPz,zz+dz/2);



% if any((TP0-DOP0-((Chl0 + C0)./Y_cp-S0*Fstable))<0) %NEW!!!
%     error('Sum of initial DOP, stably particle bound P, and P contained in Chl (both groups) a cannot be larger than TP')
% end

% if any((TP0-DOP0-((Chl0 + C0)./Y_cp-S0*Fstable))<0) %NEW!!!
%     error('Sum of initial DOP, stably particle bound P, and P contained in Chl_sed a cannot be larger than TP_sed')
% end

% if (any(FIM0<0)||any(FIM0>1))
%     error('Initial fraction of inorganic matter in sediments must be between 0 and 1')
% end

% if (any(ksw>(H_sed*(1-F_sed_sld))))
%     error('Parameter ksw is larger than the volume (thickness) of porewater')
% end  %OBS! Ideally should also be that the daily diffused porewater should not be larger
% %than the corresponding water layer volume, but this seems very unlike in practise

Tz = T0;
Cz = C0; % (mg m-3)
POCz = POC0; % (kg m-3)
Chlz = Chl0;  % (mg m-3)
DOPz = DOP0;  % (mg m-3)
DOCz = DOC0;   % (mg m-3)
DICz = DIC0;   % (mg m-3)
O2z = O20;   % (mg m-3)

NO3z = NO30;
NH4z =NH40;
SO4z = SO40;
HSz = HS0;
H2Sz = H2S0;
Fe2z = Fe20;
Ca2z = Ca20;
pHz = pH0;
CH4z = CH40;
Fe3z = Fe30;
Al3z = Al30;
SiO4z = SiO40;
SiO2z = SiO20;
diatomz = diatom0;
POPz = POP0;
Pz = (TP0-DOP0-POP0) / 2;
PPz = (TP0-DOP0-POP0) / 2; % (mg m-3) NEW!!!
Pz = Pz .* (Pz>0);
PPz = PPz .* (PPz>0);

% assume linear initial temperature profile in sediment (4 deg C at the bottom)
clear Tzy_sed
for j=1:Nz
    Tzy_sed(:,j) = interp1([0.2 10], [Tz(j) 4], [0.2:0.2:2 2.5:0.5:10])';
end

S_resusp=S_res_hypo*ones(Nz,1); %hypolimnion resuspension assumed on the first time step

rho_snow=rho_new_snow;   %initial snow density (kg m-3)
Tice=NaN;                %ice surface temperature (initial value, deg C)
XE_melt=0;               %energy flux that is left from last ice melting (initial value, W m-2)
XE_surf=0;               %energy flux from water to ice (initial value,  J m-2 per day)

%Initialisation of ice & snow variables
Hi=Ice0(1);               %total ice thickness (initial value, m)
WEQs=(rho_snow/rho_fw)*Ice0(2); %snow water equivalent  (initial value, m)
Hsi=0;                %snow ice thickness (initial value = 0 m)
HFrazil=0;              % (initial value, m) NEW!!!


if ((Hi<=0)&&(WEQs>0))
    error('Mismatch in initial ice and snow thicknesses')
end

if (Hi<=0)
    IceIndicator=0;     %IceIndicator==0 means no ice cover
else
    IceIndicator=1;
end

pp=1; %initial indexes for ice freezing/melting date arrays
qq=1;
DoF=[]; %initialize
DoM=[]; %initialize

% ============ sediment module ============
% Allocation and initial sediment profiles concentrations and reading initial concentrations for sediment from file
[sediment_concentrations, sediment_params, sediment_matrix_templates]  = sediment_init( pH, zm, In_Tz(end) );

% Passing MyLake parameters in Chemical module
mylake_params.dz = dz; mylake_params.zm = zm; mylake_params.zz = zz; mylake_params.Kz_K1 = Kz_K1; mylake_params.Kz_K1_ice = Kz_K1_ice; mylake_params.Kz_N0 = Kz_N0; mylake_params.C_shelter = C_shelter; mylake_params.lat = lat; mylake_params.lon = lon; mylake_params.alb_melt_ice = alb_melt_ice; mylake_params.alb_melt_snow = alb_melt_snow; mylake_params.PAR_sat = PAR_sat; mylake_params.f_par = f_par; mylake_params.beta_chl = beta_chl; mylake_params.lambda_i = lambda_i; mylake_params.lambda_s = lambda_s; mylake_params.F_sed_sld = F_sed_sld; mylake_params.I_scV = I_scV; mylake_params.I_scT = I_scT; mylake_params.I_scC = I_scC; mylake_params.I_scPOC = I_scPOC; mylake_params.I_scTP = I_scTP; mylake_params.I_scDOP = I_scDOP; mylake_params.I_scChl = I_scChl; mylake_params.I_scDOC = I_scDOC; mylake_params.I_scPOP = I_scPOP; mylake_params.I_scO = I_scO; mylake_params.I_scDIC = I_scDIC; mylake_params.I_scNO3 = I_scNO3; mylake_params.I_scNH4 = I_scNH4; mylake_params.I_scSO4 = I_scSO4; mylake_params.I_scFe2 = I_scFe2; mylake_params.I_scCa2 = I_scCa2; mylake_params.I_scpH = I_scpH; mylake_params.I_scCH4 = I_scCH4; mylake_params.I_scFe3 = I_scFe3; mylake_params.I_scAl3 = I_scAl3; mylake_params.I_scSiO4 = I_scSiO4; mylake_params.I_scSiO2 = I_scSiO2; mylake_params.I_scdiatom = I_scdiatom; mylake_params.swa_b0 = swa_b0; mylake_params.swa_b1 = swa_b1; mylake_params.S_res_epi = S_res_epi; mylake_params.S_res_hypo = S_res_hypo; mylake_params.H_sed = H_sed; mylake_params.Psat_L = Psat_L; mylake_params.Fmax_L = Fmax_L; mylake_params.w_s = w_s; mylake_params.w_chl = w_chl; mylake_params.Y_cp = Y_cp; mylake_params.m_twty = m_twty; mylake_params.g_twty = g_twty; mylake_params.k_twty = k_twty; mylake_params.dop_twty = dop_twty; mylake_params.P_half = P_half; mylake_params.PAR_sat_2 = PAR_sat_2; mylake_params.beta_chl_2 = beta_chl_2; mylake_params.w_chl_2 = w_chl_2; mylake_params.m_twty_2 = m_twty_2; mylake_params.g_twty_2 = g_twty_2; mylake_params.P_half_2 = P_half_2; mylake_params.oc_DOC = oc_DOC; mylake_params.qy_DOC = qy_DOC; mylake_params.k_BOD = k_BOD; mylake_params.k_SOD = k_SOD; mylake_params.theta_bod = theta_bod; mylake_params.theta_bod_ice = theta_bod_ice; mylake_params.theta_sod = theta_sod; mylake_params.theta_sod_ice = theta_sod_ice; mylake_params.BOD_temp_switch = BOD_temp_switch; mylake_params.pH = pH; mylake_params.Q10 = Q10; mylake_params.wc_factor = wc_factor; mylake_params.T_ref = T_ref; mylake_params.theta_m = theta_m;  mylake_params.floculation_switch = floculation_switch; mylake_params.rate_estimator_switch = rate_estimator_switch; mylake_params.Az = Az; mylake_params.Vz = Vz; mylake_params.dt = dt;


% >>>>>> Start of the time loop >>>>>>
Resuspension_counter=zeros(Nz,1); %kg
Sedimentation_counter=zeros(Nz,1); %kg
SS_decr=0; %kg

for i = 1:length(tt)

    % Surface heat fluxes (W m-2), wind stress (N m-2) & daylight fraction (-), based on Air-Sea Toolbox
    [Qsw,Qlw,Qsl,tau,DayFrac,DayFracHeating] = heatflux_v12(tt(i),Wt(i,1),Wt(i,2),Wt(i,3),Wt(i,4),Wt(i,5),Wt(i,6),Tz(1), ...
        lat,lon,WEQs,Hi,alb_melt_ice,alb_melt_snow,albedot1);     %Qlw and Qsl are functions of Tz(1)

    % Calculate total mean PAR and non-PAR light extinction coefficient in water (background + due to Chl a)
    lambdaz_wtot_avg=zeros(Nz,1);
    lambdaz_NP_wtot_avg=zeros(Nz,1);

    %NEW!!! below additional term for chlorophyll group 2
    if (selfshading_switch==1)
        lambdaz_wtot=swa_b1 * ones(Nz,1) + beta_chl*Chlz + beta_chl_2*Cz; %at layer z
        lambdaz_NP_wtot=swa_b0 * ones(Nz,1) + beta_chl*Chlz + beta_chl_2*Cz; %at layer z
        for j=1:Nz
            lambdaz_wtot_avg(j)=mean(swa_b1 * ones(j,1) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
            lambdaz_NP_wtot_avg(j)=mean(swa_b0 * ones(j,1) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
        end
    else %constant with depth
        lambdaz_wtot=swa_b1 * ones(Nz,1);
        lambdaz_wtot_avg=swa_b1 * ones(Nz,1);
        lambdaz_NP_wtot=swa_b0 * ones(Nz,1);
        lambdaz_NP_wtot_avg=swa_b0 * ones(Nz,1);
    end %if selfshading...


    if(IceIndicator==0)
        IceSnowAttCoeff=1; %no extra light attenuation due to snow and ice
    else    %extra light attenuation due to ice and snow
        IceSnowAttCoeff=exp(-lambda_i * Hi) * exp(-lambda_s * (rho_fw/rho_snow)*WEQs);
    end

    Tprof_prev=Tz; %temperature profile at previous time step (for convection_v2.m)

    rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);  % Density (kg/m3)

    % Sediment vertical heat flux, Q_sed
    % (averaged over the whole top area of the layer, although actually coming only from the "sides")
    if (sediment_heatflux_switch==1)
        % update top sediment temperatures
        dz_sf = 0.2; %fixed distance between the two topmost sediment layers (m)
        Tzy_sed(1,:) = Tz';
        Tzy_sed_upd = sedimentheat_v11(Tzy_sed, K_sed, dt);
        Tzy_sed=Tzy_sed_upd;
        Qz_sed=K_sed*rho_sed*cp_sed*(1/dz_sf)*(-diff([Az; 0])./Az) .* (Tzy_sed(2,:)'-Tzy_sed(1,:)'); %(J day-1 m-2)
        %positive heat flux => from sediment to water
    else
        Qz_sed = zeros(Nz,1);
    end

    Cw = 4.18e+6;	% Volumetric heat capacity of water (J K-1 m-3)

    %Heat sources/sinks:
    %Total attenuation coefficient profile, two-band extinction, PAR & non-PAR
    Par_Attn=exp([0; -lambdaz_wtot_avg] .* [zz; zz(end)+dz]);
    NonPar_Attn=exp([0; -lambdaz_NP_wtot_avg] .* [zz; zz(end)+dz]);

    Attn_z=(-f_par * diff([1; ([Az(2:end);0]./Az).*Par_Attn(2:end)]) + ...
        (-(1-f_par)) * diff([1; ([Az(2:end);0]./Az).*NonPar_Attn(2:end)])); %NEW (corrected 210807)

    if(IceIndicator==0)
        % 1) Vertical heating profile for open water periods (during daytime heating)
        Qz = (Qsw + XE_melt) * Attn_z; %(W m-2)
        Qz(1) = Qz(1) + DayFracHeating*(Qlw + Qsl); %surface layer heating
        XE_melt=0; %Reset
        dT = Az .* ((60*60*24*dt) * Qz + DayFracHeating*Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (daytime heating, ice melt, sediment);

        % === Frazil ice melting, NEW!!! === %
        postemp=find(dT>0);
        if (isempty(postemp)==0)
            RelT=dT(postemp)./sum(dT(postemp));
            HFrazilnew=max(0, HFrazil - sum(dT(postemp))*1/((Az(1)*rho_ice*L_ice)/(Cw * Vz(1)))); %
            sumdTnew = max(0, sum(dT(postemp))-(HFrazil*Az(1)*rho_ice*L_ice)/(Cw * Vz(1)));
            dT(postemp)=RelT.*sumdTnew;
            HFrazil=HFrazilnew;
        end
        % === === ===
    else
        % Vertical heating profile for ice-covered periods (both day- and nighttime)
        Qz = Qsw * IceSnowAttCoeff * Attn_z; %(W/m2)
        dT = Az .* ((60*60*24*dt) * Qz + Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (solar rad., sediment);
    end

    Tz = Tz + dT;        %Temperature change after daytime surface heatfluxes (or whole day in ice covered period)

    % Convective mixing adjustment (mix successive layers until stable density profile)
    % and
    % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
    [Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz] = convection_v2(Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1);
    Tprof_prev=Tz; %NEW!!! Update Tprof_prev

    if(IceIndicator==0)
        % 2) Vertical heating profile for open water periods (during nighttime heating)
        [Qsw,Qlw_2,Qsl_2,tau,DayFrac,DayFracHeating] = heatflux_v12(tt(i),Wt(i,1),Wt(i,2),Wt(i,3),Wt(i,4),Wt(i,5),Wt(i,6),Tz(1), ...
            lat,lon,WEQs,Hi,alb_melt_ice,alb_melt_snow,albedot1); %Qlw and Qsl are functions of Tz(1)
        Qz(1) = (1-DayFracHeating)*(Qlw_2 + Qsl_2); %surface layer heating
        Qz(2:end)=0; %No other heating below surface layer
        dT = Az .* ((60*60*24*dt) * Qz + (1-DayFracHeating)*Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (longwave & turbulent fluxes);

        % === NEW!!! frazil ice melting === %
        postemp=find(dT>0);
        if (isempty(postemp)==0)
            %disp(['NOTE: positive night heat flux at T=' num2str(Tz(postemp),2)]) %NEW
            RelT=dT(postemp)./sum(dT(postemp));
            HFrazilnew=max(0, HFrazil - sum(dT(postemp))*1/((Az(1)*rho_ice*L_ice)/(Cw * Vz(1)))); %
            sumdTnew = max(0, sum(dT(postemp))-(HFrazil*Az(1)*rho_ice*L_ice)/(Cw * Vz(1)));
            dT(postemp)=RelT.*sumdTnew;
            HFrazil=HFrazilnew;
        end
        % === === ===

        Tz = Tz + dT;         %Temperature change after nighttime surface heatfluxes

        % Convective mixing adjustment (mix successive layers until stable density profile)
        % and
        % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
        [Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz] = convection_v2(Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1);

        Qlw = DayFracHeating*Qlw + (1-DayFracHeating)*Qlw_2; %total amounts, only for output purposes
        Qsl = DayFracHeating*Qsl + (1-DayFracHeating)*Qsl_2; %total amounts, only for output purposes
    end

    % Vertical turbulent diffusion
    g   = 9.81;							% Gravity acceleration (m s-2)
    rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);  % Water density (kg m-3)
    % Note: in equations of rho it is assumed that every supercooled degree lowers density by
    % 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")

    N2  = g * (diff(log(rho)) ./ diff(zz));	% Brunt-Vaisala frequency (s-2) for level (zz+1)
    if (IceIndicator==0)
        Kz  = Kz_K1 * max(Kz_N0, N2).^(-Kz_b1);	% Vertical diffusion coeff. in ice free season (m2 day-1)
        % for level (zz+1)
    else
        Kz  = Kz_K1_ice * max(Kz_N0, N2).^(-Kz_b1_ice); % Vertical diffusion coeff. under ice cover (m2 day-1)
        % for level (zz+1)
    end

    Fi = tridiag_DIF_v11([NaN; Kz],Vz,Az,dz,dt); %Tridiagonal matrix for general diffusion

    Tz = Fi \ (Tz);        %Solving new temperature profile (diffusion, sources/sinks already added to Tz above)

    % Convective mixing adjustment (mix successive layers until stable density profile)
    % (don't allow temperature jumps over temperature of maximum density, no summer/autumn turnover here!)
    [Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz] = convection_v2(Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,0);

    %% Deposition
    if deposition_switch == 1
        Pz(1) = Pz(1) + (Deposition(i,5) ./ sum(Vz(1))) ; % qty added in mg to the top layer z
    end

    % NEW!!! === Code rearranging
    % Calculate again the total mean PAR light extinction coefficient in water (background + due to Chl a)
    lambdaz_wtot_avg=zeros(Nz,1);

    %NEW!!! below additional term for chlorophyll group 2
    if (selfshading_switch==1)
        lambdaz_wtot=swa_b1 * ones(Nz,1) + beta_chl*Chlz + beta_chl_2*Cz; %at layer z.
        for j=1:Nz
            lambdaz_wtot_avg(j)=mean(swa_b1 * ones(j,1) + beta_chl*Chlz(1:j) + beta_chl_2*Cz(1:j)); %average down to layer z
        end
    else %constant with depth
        lambdaz_wtot=swa_b1 * ones(Nz,1);
        lambdaz_wtot_avg=swa_b1 * ones(Nz,1);
    end %if selfshading...

    %Photosynthetically Active Radiation (for chlorophyll group 1)
    H_sw_z=NaN*zeros(Nz,1);

    % ===== NEW!!! bug (when Dayfrac==0) fixed 071107
    if ((IceIndicator==0)&&(DayFrac>0))
        PAR_z=((3/2) / (e_par * DayFrac)) * f_par * Qsw  * exp(-lambdaz_wtot_avg .* zz);
        %Irradiance at noon (mol m-2 s-1) at levels zz
    elseif ((IceIndicator==1)&&(DayFrac>0))    %extra light attenuation due to ice and snow
        PAR_z=((3/2) / (e_par * DayFrac)) * IceSnowAttCoeff * f_par *...
            Qsw  * exp(-lambdaz_wtot_avg .* zz);
    else PAR_z=zeros(Nz,1); %DayFrac==0, polar night
    end
    % =====
    U_sw_z=PAR_z./PAR_sat; %scaled irradiance at levels zz
    inx_u=find(U_sw_z<=1); %undersaturated
    inx_s=find(U_sw_z>1);  %saturated

    H_sw_z(inx_u)=(2/3)*U_sw_z(inx_u);  %undersaturated

    dum_a=sqrt(U_sw_z);
    dum_b=sqrt(U_sw_z-1);
    H_sw_z(inx_s)=(2/3)*U_sw_z(inx_s) + log((dum_a(inx_s) + dum_b(inx_s))./(dum_a(inx_s) ...  %saturated
        - dum_b(inx_s))) - (2/3)*(U_sw_z(inx_s)+2).*(dum_b(inx_s)./dum_a(inx_s));


    %Photosynthetically Active Radiation (for chlorophyll group 2) NEW!!!
    H_sw_z_2=NaN*zeros(Nz,1);

    U_sw_z_2=PAR_z./PAR_sat_2; %scaled irradiance at levels zz
    inx_u=find(U_sw_z_2<=1); %undersaturated
    inx_s=find(U_sw_z_2>1);  %saturated

    H_sw_z_2(inx_u)=(2/3)*U_sw_z_2(inx_u);  %undersaturated

    dum_a=sqrt(U_sw_z_2);
    dum_b=sqrt(U_sw_z_2-1);
    H_sw_z_2(inx_s)=(2/3)*U_sw_z_2(inx_s) + log((dum_a(inx_s) + dum_b(inx_s))./(dum_a(inx_s) ...  %saturated
        - dum_b(inx_s))) - (2/3)*(U_sw_z_2(inx_s)+2).*(dum_b(inx_s)./dum_a(inx_s));


    % NOTE: All reactions are moved in "rates" method below on 26.03.2017
    % Corrected on 26.06.2017: diffusion and advection for all species in WC

    Fi_ad_w_s = tridiag_HAD_v11([NaN; Kz],w_s,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion
    Fi_ad_w_chl = tridiag_HAD_v11([NaN; Kz],w_chl,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion
    Fi_ad_w_chl_2 = tridiag_HAD_v11([NaN; Kz],w_chl_2,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion

    Cz = Fi_ad_w_chl_2 \ (Cz);  %Solving new phytoplankton profile (advection + diffusion) (always larger than background level)
    POCz = Fi_ad_w_s \ (POCz);           %Solving new suspended solids profile (advection + diffusion)
    Pz = Fi \ (Pz); %Solving new dissolved inorganic P profile (diffusion)
    Chlz = Fi_ad_w_chl \ (Chlz);  %Solving new phytoplankton profile (advection + diffusion) (always larger than background level)
    PPz = Fi_ad_w_s \ (PPz);     %Solving new suspended particulate inorganic P profile (advection + diffusion)
    DOPz = Fi \ (DOPz);
    % DOCz  below
    % DICz  below
    % O2z  below
    NO3z = Fi \ NO3z;
    NH4z = Fi \ NH4z;
    SO4z = Fi \ SO4z;
    HSz = Fi \ HSz;
    H2Sz = Fi \ H2Sz;
    Fe2z = Fi \ Fe2z;
    Ca2z = Fi \ Ca2z;
    Fe3z = Fi_ad_w_s \ (Fe3z);
    Al3z = Fi_ad_w_s \ (Al3z);
    SiO4z = Fi \ (SiO4z);
    SiO2z = Fi_ad_w_s \ (SiO2z);
    diatomz = Fi_ad_w_s \ (diatomz);
    POPz = Fi_ad_w_s \ (POPz);


    %Dissolved organic carbon
    % - current version
    Kd_old=0;
    Theeta=0;
    Currdate=datevec(tt(i)); %Date
    %Date=0;
    Date=Currdate(1,2); %Month number
    if (photobleaching==1) %Fokema
        %[DOCz,Kd_new] = fokema(DOCz,Kd_old,Qsw,Tz,Theeta,Date,zz);
        DOCz1 = 0.0775.*DOCz; %Subpools
        DOCz2 = 0.1486.*DOCz;
        DOCz3 = 0.7739.*DOCz;
        [DOCz1_new,DOCz2_new,DOCz3_new,DOC1frac,DOC2frac,DOC3frac,Kd_new,Daily_BB1,Daily_BB2,Daily_BB3,Daily_PB] = fokema_new(DOCz1,DOCz2,DOCz3,Kd_old,Qsw,Tz,Theeta,Date,zz);
        DOCz = DOCz1_new + DOCz2_new + DOCz3_new; %Total DOC
        DOCz = Fi \ DOCz; %Solving new dissolved organic C profile (diffusion)
        %DOCz1_new = Fi \ DOCz1_new; %Solving new dissolved organic C profile (diffusion)
        %DOCz2_new = Fi \ DOCz2_new; %Solving new dissolved organic C profile (diffusion)
        %DOCz3_new = Fi \ DOCz3_new; %Solving new dissolved organic C profile (diffusion)

    else %TSA model
        dDOC = -oc_DOC*qy_DOC*f_par*(1/e_par)*(60*60*24*dt)*Qsw*Attn_z; %photochemical degradation
        %[m2/mg_doc]*[mg_doc/mol_qnt]*[-]*[mol_qnt/J]*[s/day]*[J/s/m2]*[-] = [1/day]
        DOCz = Fi \ (DOCz + dDOC.*DOCz); %Solving new dissolved organic C profile (diffusion)
    end


    %Oxygen surface flux
    if(IceIndicator==0)
        [O2z(1),O2flux,O2_eq,K0_O2] = oxygenflux(O2z(1),C_shelter^(1/3)*Wt(i,6),Wt(i,5),Tz(1),dz);
    else
        O2flux = 0;
    end

    O2z = Fi \ O2z; %Solving new dissolved oxygen profile (diffusion)
    DOCz = Fi \ DOCz;

    %Dissolved inorganic carbon
    %DIC partitioning in water
    [CO2z,CO2frac] = carbonequilibrium(DICz,Tz,pH);


    % CO2 production by degraded DOC
    % CO2z = max(0,CO2z + 1.375.*(-O2_diff));
    DICz = CO2z./CO2frac;
    %TC = Tz(1); %For monitoring only

    %Carbon dioxide surface flux
    if(IceIndicator==0)
        [CO2z(1),surfflux,CO2_eq,K0,CO2_ppm] = carbondioxideflux(CO2z(1),C_shelter^(1/3)*Wt(i,6),Wt(i,5),Tz(1),dz,tt(i));
        DICz(1) = CO2z(1)/CO2frac(1);
    else
        surfflux=0;
    end

    DICz = Fi \ DICz; %Solving new DIC profile (diffusion)

    %==================================
    %Dissolved inorganic carbon - version II
    %
    %     pH = 7*ones(1,length(zz));
    %     TC = Tz(1); %For monitoring only
    %     [CO2z,CO2frac,C_acid,C_basic] = carbonequilibrium(DICz,Tz,pH);
    %
    %     C_acid = C_acid + 1.375.*(-O2_erotus);
    %     DICz = C_acid+C_basic;
    %     [CO2z,CO2frac,C_acid,C_basic] = carbonequilibrium(DICz,Tz,pH);
    %     if(IceIndicator==0)
    %         [CO2z(1),surfflux,CO2_eq,K0,CO2_ppm] = carbondioxideflux(CO2z(1),Wt(i,6),Wt(i,5),Tz(1),dz,tt(i));
    %         DICz(1) = CO2z(1)/CO2frac(1);
    %     else
    %         surfflux=0;
    %     end
    %     DICz = Fi \ DICz;
    %==================================

    %Sediment-water exchange (DOP source neglected)
    %-porewater to water

    if resuspension_enabled == 0
        ksw = 0;
    end


    if (river_inflow_switch==1)
        Iflw = I_scV * Inflw(i,1); % (scaled) inflow rate
        Iflw_T = I_scT + Inflw(i,2); %(adjusted) inflow temperature
        if (Iflw_T<Tf) %negative temperatures changed to Tf
            Iflw_T=Tf;
        end
        Iflw_C = I_scC * Inflw(i,3); %(scaled) inflow C concentration
        Iflw_POC = I_scPOC * Inflw(i,4); %(scaled) inflow POC concentration
        Iflw_TP = I_scTP * Inflw(i,5); %(scaled) inflow TP concentration (incl. DOP & Chla)
        Iflw_DOP = I_scDOP * Inflw(i,6); %(scaled) inflow DOP concentration
        Iflw_Chl = I_scChl * Inflw(i,7); %(scaled) inflow Chl a concentration
        Iflw_DOC = I_scDOC * Inflw(i,8); %(scaled) inflow DOC concentration
        Iflw_DIC = I_scDIC * Inflw(i,9); %(scaled) inflow DIC concentration
        Iflw_O2 = I_scO * Inflw(i,10); %(scaled) inflow O2 concentration

        % inflow HS and H2S are neglected

        % Where inflow POPC???

        Iflw_NO3 = I_scNO3 * Inflw(i,11);
        Iflw_NH4 = I_scNH4 * Inflw(i,12);
        Iflw_SO4 = I_scSO4 * Inflw(i,13);
        Iflw_Fe2 = I_scFe2 * Inflw(i,14);
        Iflw_Ca2 = I_scCa2 * Inflw(i,15);
        Iflw_pH = I_scpH * Inflw(i,16);
        Iflw_CH4 = I_scCH4 * Inflw(i,17);
        Iflw_Fe3 = I_scFe3 * Inflw(i,18);
        Iflw_Al3 = I_scAl3 * Inflw(i,19);
        Iflw_SiO4 = I_scSiO4 * Inflw(i,20);
        Iflw_SiO2 = I_scSiO2 * Inflw(i,21);
        Iflw_diatom = I_scdiatom * Inflw(i,22);
        Iflw_POP = I_scPOP * Inflw(i,23);
        Iflw_Pz = (Iflw_TP - Iflw_DOP - Iflw_POP)/2;
        Iflw_Pz = Iflw_Pz .* (Iflw_Pz > 0);
        Iflw_PP = (Iflw_TP - Iflw_DOP - Iflw_POP)/2; %;
        Iflw_PP = Iflw_PP .* (Iflw_PP > 0);
        % Iflw_PP = (Iflw_TP - Iflw_DOP - Iflw_Chl - Iflw_Chl)/2;
        % Iflw_Pz = (Iflw_TP - Iflw_DOP - Iflw_Chl - Iflw_Chl)/2;


        % %Added suspended solids correction: minimum allowed P bioavailability factor is 0.1
        % if any((1-(Iflw_DOP+(Iflw_Chl+Iflw_C)./Y_cp)./Iflw_TP-(Iflw_S*Fstable)./Iflw_TP)<0.1); % NEW!!!!
        %     Iflw_S_dum = (1 - (Iflw_DOP+(Iflw_Chl+Iflw_C)./Y_cp)./Iflw_TP - 0.1).*(Iflw_TP./Fstable); %NEW!!!
        %     SS_decr=SS_decr+(Iflw_S-Iflw_S_dum)*Iflw;
        %     Iflw_S=Iflw_S_dum;
        % end

        % if any((Iflw_TP-Iflw_DOP-(Iflw_Chl+Iflw_C)./Y_cp-Iflw_S*Fstable)<0)  %NEW!!!
        %     error('Sum of DOP, inactive PP, and P contained in Chl a (both groups) in inflow cannot be larger than TP')
        % end


        if(Iflw>0)
            if (isnan(Iflw_T))
                lvlD=0;
                Iflw_T=Tz(1);
            else
                rho = polyval(ies80,max(0,Tz(:)))+min(Tz(:),0);	% Density (kg/m3)
                rho_Iflw=polyval(ies80,max(0,Iflw_T))+min(Iflw_T,0);
                lvlG=find(rho>=rho_Iflw);
                if (isempty(lvlG))
                    lvlG=length(rho);
                end
                lvlD=zz(lvlG(1)); %level zz above which inflow is put
            end %if isnan...


            %Changes in properties due to inflow
            Tz=IOflow_v11(dz, zz, Vz, Tz, lvlD, Iflw, Iflw_T); %Temperature
            POCz=IOflow_v11(dz, zz, Vz, POCz, lvlD, Iflw, Iflw_POC); % POC
            DOPz=IOflow_v11(dz, zz, Vz, DOPz, lvlD, Iflw, Iflw_DOP); %Particulate organic P

            % TODO: this needs to be moved in the reaction module too.
            % TIPz=Pz + PPz; % Total inorg. phosphorus (excl. Chla and DOP) in the water column (mg m-3)

             %Total inorg. phosphorus (excl. Chla and DOP)
            % TIPz=IOflow_v11(dz, zz, Vz, TIPz, lvlD, Iflw, Iflw_TP-((Iflw_Chl+Iflw_C)./Y_cp)-Iflw_DOP); %NEW!!!


            %== P-partitioning in water==
            % [Pz, trash]=Ppart(POCz./rho_sed,TIPz,Psat_L,Fmax_L,rho_sed,Fstable);
            % PPz=TIPz-Pz;

            Cz=IOflow_v11(dz, zz, Vz, Cz, lvlD, Iflw, Iflw_C); %Chlorophyll (group 2)
            Chlz=IOflow_v11(dz, zz, Vz, Chlz, lvlD, Iflw, Iflw_Chl); %Chlorophyll (group 1)
            DOCz=IOflow_v11(dz, zz, Vz, DOCz, lvlD, Iflw, Iflw_DOC); %DOC
            DICz=IOflow_v11(dz, zz, Vz, DICz, lvlD, Iflw, Iflw_DIC); %DIC
            O2z=IOflow_v11(dz, zz, Vz, O2z, lvlD, Iflw, Iflw_O2); %O2
            NO3z=IOflow_v11(dz, zz, Vz, NO3z, lvlD, Iflw, Iflw_NO3); %NO3
            NH4z=IOflow_v11(dz, zz, Vz, NH4z, lvlD, Iflw, Iflw_NH4); %NH4
            SO4z=IOflow_v11(dz, zz, Vz, SO4z, lvlD, Iflw, Iflw_SO4); %SO4
            Fe2z=IOflow_v11(dz, zz, Vz, Fe2z, lvlD, Iflw, Iflw_Fe2); %Fe2
            Ca2z=IOflow_v11(dz, zz, Vz, Ca2z, lvlD, Iflw, Iflw_Ca2); %Ca2
            pHz=IOflow_v11(dz, zz, Vz, pHz, lvlD, Iflw, Iflw_pH); %pH
            CH4z=IOflow_v11(dz, zz, Vz, CH4z, lvlD, Iflw, Iflw_CH4); %CH4
            Fe3z=IOflow_v11(dz, zz, Vz, Fe3z, lvlD, Iflw, Iflw_Fe3); %Fe3
            Al3z=IOflow_v11(dz, zz, Vz, Al3z, lvlD, Iflw, Iflw_Al3); %Al3
            SiO4z=IOflow_v11(dz, zz, Vz, SiO4z, lvlD, Iflw, Iflw_SiO4); %SiO4
            SiO2z=IOflow_v11(dz, zz, Vz, SiO2z, lvlD, Iflw, Iflw_SiO2); %SiO2
            diatomz=IOflow_v11(dz, zz, Vz, diatomz, lvlD, Iflw, Iflw_diatom); %diatom
            POPz=IOflow_v11(dz, zz, Vz, POPz, lvlD, Iflw, Iflw_POP); %POP
            PPz=IOflow_v11(dz, zz, Vz, PPz, lvlD, Iflw, Iflw_PP); %PP
            Pz=IOflow_v11(dz, zz, Vz, Pz, lvlD, Iflw, Iflw_Pz); %Pz
        else
            lvlD=NaN;
        end %if(Iflw>0)

    else
        Iflw=0; % (scaled) inflow rate
        Iflw_T = NaN; %(adjusted) inflow temperature
        Iflw_C = NaN; %(scaled) inflow C concentration
        Iflw_POC = NaN; %(scaled) inflow S concentration
        Iflw_TP = NaN; %(scaled) inflow TP concentration (incl. DOP & Chla)
        Iflw_DOP = NaN; %(scaled) inflow DOP concentration
        Iflw_Chl = NaN; %(scaled) inflow Chl a concentration
        Iflw_DOC = NaN; %(scaled) inflow DOC concentration
        Iflw_DIC = NaN; %(scaled) inflow DIC concentration
        Iflw_O2 = NaN; %(scaled) inflow concentration
        Iflw_NO3 = NaN; %(scaled) inflow concentration
        Iflw_NH4 = NaN; %(scaled) inflow concentration
        Iflw_SO4 = NaN; %(scaled) inflow  concentration
        Iflw_Fe2 = NaN; %(scaled) inflow concentration
        Iflw_Ca2 = NaN; %(scaled) inflow concentration
        Iflw_pH = NaN; %(scaled) inflow concentration
        Iflw_CH4 = NaN; %(scaled) inflow concentration
        Iflw_Fe3 = NaN; %(scaled) inflow concentration
        Iflw_Al3 = NaN; %(scaled) inflow concentration
        Iflw_SiO4 = NaN; %(scaled) inflow concentration
        Iflw_SiO2 = NaN; %(scaled) inflow concentration
        Iflw_diatom = NaN; %(scaled) inflow concentration
        lvlD=NaN;
    end  %if (river_inflow_switch==1)

    % Convective mixing adjustment (mix successive layers until stable density profile,  no summer/autumn turnover here!)

    [Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz] = convection_v2(Tz,Cz,POCz,Pz,Chlz,PPz,DOPz,DOCz,DICz,O2z,NO3z,NH4z,SO4z,HSz,H2Sz,Fe2z,Ca2z,pHz,CH4z,Fe3z,Al3z,SiO4z,SiO2z,diatomz,POPz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,0);

    if (IceIndicator==0)

        TKE=C_shelter*Az(1)*sqrt(tau^3/rho(1))*(24*60*60*dt); %Turbulent kinetic energy (J day-1) over the whole lake

        %Wind mixing
        WmixIndicator=1;
        Bef_wind=sum(diff(rho)==0); %just a watch variable
        while (WmixIndicator==1)
            d_rho=diff(rho);
            inx=find(d_rho>0);
            if (isempty(inx)==0); %if water column not already fully mixed
                zb=inx(1);
                MLD=dz*zb; %mixed layer depth
                dD=d_rho(zb); %density difference
                Zg=sum( Az(1:zb+1) .* zz(1:zb+1) ) / sum(Az(1:zb+1)); %Depth of center of mass of mixed layer
                V_weight=Vz(zb+1)*sum(Vz(1:zb))/(Vz(zb+1)+sum(Vz(1:zb)));
                POE=(dD*g*V_weight*(MLD + dz/2 - Zg));
                KP_ratio=TKE/POE;
                if (KP_ratio>=1)

                    Tmix=sum( Vz(1:zb+1).*Tz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Tz(1:zb+1)=Tmix;

                    Cmix=sum( Vz(1:zb+1).*Cz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Cz(1:zb+1)=Cmix;

                    Smix=sum( Vz(1:zb+1).*POCz(1:zb+1) ) / sum(Vz(1:zb+1));
                    POCz(1:zb+1)=Smix;

                    Pmix=sum( Vz(1:zb+1).*Pz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Pz(1:zb+1)=Pmix;

                    Chlmix=sum( Vz(1:zb+1).*Chlz(1:zb+1) ) / sum(Vz(1:zb+1));
                    Chlz(1:zb+1)=Chlmix;

                    PPmix=sum( Vz(1:zb+1).*PPz(1:zb+1) ) / sum(Vz(1:zb+1));
                    PPz(1:zb+1)=PPmix;

                    DOPmix=sum( Vz(1:zb+1).*DOPz(1:zb+1) ) / sum(Vz(1:zb+1));
                    DOPz(1:zb+1)=DOPmix;

                    DOCmix=sum( Vz(1:zb+1).*DOCz(1:zb+1) ) / sum(Vz(1:zb+1));
                    DOCz(1:zb+1)=DOCmix;

                    DICmix=sum( Vz(1:zb+1).*DICz(1:zb+1) ) / sum(Vz(1:zb+1));
                    DICz(1:zb+1)=DICmix;

                    O2mix=sum( Vz(1:zb+1).*O2z(1:zb+1) ) / sum(Vz(1:zb+1));
                    O2z(1:zb+1)=O2mix;

                    NO3mix=sum( Vz(1:zb+1).*NO3z(1:zb+1) ) / sum(Vz(1:zb+1));
                    NO3z(1:zb+1)=NO3mix;

                    NH4mix=sum( Vz(1:zb+1).*NH4z(1:zb+1) ) / sum(Vz(1:zb+1));
                    NH4z(1:zb+1)=NH4mix;

                    SO4mix=sum( Vz(1:zb+1).*SO4z(1:zb+1) ) / sum(Vz(1:zb+1));
                    SO4z(1:zb+1)=SO4mix;

                    HSmix=sum( Vz(1:zb+1).*HSz(1:zb+1) ) / sum(Vz(1:zb+1));
                    HSz(1:zb+1)=HSmix;

                    H2Smix=sum( Vz(1:zb+1).*H2Sz(1:zb+1) ) / sum(Vz(1:zb+1));
                    H2Sz(1:zb+1)=H2Smix;

                    Fe2mix=sum( Vz(1:zb+1).*Fe2z(1:zb+1) ) / sum(Vz(1:zb+1));
                    Fe2z(1:zb+1)=Fe2mix;

                    Ca2mix=sum( Vz(1:zb+1).*Ca2z(1:zb+1) ) / sum(Vz(1:zb+1));
                    Ca2z(1:zb+1)=Ca2mix;

                    pHmix=sum( Vz(1:zb+1).*pHz(1:zb+1) ) / sum(Vz(1:zb+1));
                    pHz(1:zb+1)=pHmix;

                    CH4mix=sum( Vz(1:zb+1).*CH4z(1:zb+1) ) / sum(Vz(1:zb+1));
                    CH4z(1:zb+1)=CH4mix;

                    Fe3mix=sum( Vz(1:zb+1).*Fe3z(1:zb+1) ) / sum(Vz(1:zb+1));
                    Fe3z(1:zb+1)=Fe3mix;

                    Al3mix=sum( Vz(1:zb+1).*Al3z(1:zb+1) ) / sum(Vz(1:zb+1));
                    Al3z(1:zb+1)=Al3mix;

                    SiO4mix=sum( Vz(1:zb+1).*SiO4z(1:zb+1) ) / sum(Vz(1:zb+1));
                    SiO4z(1:zb+1)=SiO4mix;

                    SiO2mix=sum( Vz(1:zb+1).*SiO2z(1:zb+1) ) / sum(Vz(1:zb+1));
                    SiO2z(1:zb+1)=SiO2mix;

                    diatommix=sum( Vz(1:zb+1).*diatomz(1:zb+1) ) / sum(Vz(1:zb+1));
                    diatomz(1:zb+1)=diatommix;

                    POPmix=sum( Vz(1:zb+1).*POPz(1:zb+1) ) / sum(Vz(1:zb+1));
                    POPz(1:zb+1)=POPmix;

                    rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
                    TKE=TKE-POE;
                else %if KP_ratio < 1, then mix with the remaining TKE part of the underlying layer
                    Tmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Tz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Tz(1:zb)=Tmix;
                    Tz(zb+1)=KP_ratio*Tmix + (1-KP_ratio)*Tz(zb+1);

                    Cmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Cz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Cz(1:zb)=Cmix;
                    Cz(zb+1)=KP_ratio*Cmix + (1-KP_ratio)*Cz(zb+1);

                    Smix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*POCz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    POCz(1:zb)=Smix;
                    POCz(zb+1)=KP_ratio*Smix + (1-KP_ratio)*POCz(zb+1);

                    Pmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Pz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Pz(1:zb)=Pmix;
                    Pz(zb+1)=KP_ratio*Pmix + (1-KP_ratio)*Pz(zb+1);

                    Chlmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Chlz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Chlz(1:zb)=Chlmix;
                    Chlz(zb+1)=KP_ratio*Chlmix + (1-KP_ratio)*Chlz(zb+1);

                    PPmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*PPz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    PPz(1:zb)=PPmix;
                    PPz(zb+1)=KP_ratio*PPmix + (1-KP_ratio)*PPz(zb+1);

                    DOPmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOPz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    DOPz(1:zb)=DOPmix;
                    DOPz(zb+1)=KP_ratio*DOPmix + (1-KP_ratio)*DOPz(zb+1);

                    DOCmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOCz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    DOCz(1:zb)=DOCmix;
                    DOCz(zb+1)=KP_ratio*DOCmix + (1-KP_ratio)*DOCz(zb+1);

                    DICmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DICz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    DICz(1:zb)=DICmix;
                    DICz(zb+1)=KP_ratio*DICmix + (1-KP_ratio)*DICz(zb+1);

                    O2mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*O2z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    O2z(1:zb)=O2mix;
                    O2z(zb+1)=KP_ratio*O2mix + (1-KP_ratio)*O2z(zb+1);

                    NO3mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*NO3z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    NO3z(1:zb)=NO3mix;
                    NO3z(zb+1)=KP_ratio*NO3mix + (1-KP_ratio)*NO3z(zb+1);

                    NH4mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*NH4z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    NH4z(1:zb)=NH4mix;
                    NH4z(zb+1)=KP_ratio*NH4mix + (1-KP_ratio)*NH4z(zb+1);

                    SO4mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*SO4z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    SO4z(1:zb)=SO4mix;
                    SO4z(zb+1)=KP_ratio*SO4mix + (1-KP_ratio)*SO4z(zb+1);

                    HSmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*HSz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    HSz(1:zb)=HSmix;
                    HSz(zb+1)=KP_ratio*HSmix + (1-KP_ratio)*HSz(zb+1);

                    H2Smix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*H2Sz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    H2Sz(1:zb)=H2Smix;
                    H2Sz(zb+1)=KP_ratio*H2Smix + (1-KP_ratio)*H2Sz(zb+1);

                    Fe2mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Fe2z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Fe2z(1:zb)=Fe2mix;
                    Fe2z(zb+1)=KP_ratio*Fe2mix + (1-KP_ratio)*Fe2z(zb+1);

                    Ca2mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Ca2z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Ca2z(1:zb)=Ca2mix;
                    Ca2z(zb+1)=KP_ratio*Ca2mix + (1-KP_ratio)*Ca2z(zb+1);

                    pHmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*pHz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    pHz(1:zb)=pHmix;
                    pHz(zb+1)=KP_ratio*pHmix + (1-KP_ratio)*pHz(zb+1);

                    CH4mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*CH4z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    CH4z(1:zb)=CH4mix;
                    CH4z(zb+1)=KP_ratio*CH4mix + (1-KP_ratio)*CH4z(zb+1);

                    Fe3mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Fe3z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Fe3z(1:zb)=Fe3mix;
                    Fe3z(zb+1)=KP_ratio*Fe3mix + (1-KP_ratio)*Fe3z(zb+1);

                    Al3mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Al3z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    Al3z(1:zb)=Al3mix;
                    Al3z(zb+1)=KP_ratio*Al3mix + (1-KP_ratio)*Al3z(zb+1);

                    SiO4mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*SiO4z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    SiO4z(1:zb)=SiO4mix;
                    SiO4z(zb+1)=KP_ratio*SiO4mix + (1-KP_ratio)*SiO4z(zb+1);

                    SiO2mix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*SiO2z(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    SiO2z(1:zb)=SiO2mix;
                    SiO2z(zb+1)=KP_ratio*SiO2mix + (1-KP_ratio)*SiO2z(zb+1);

                    diatommix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*diatomz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    diatomz(1:zb)=diatommix;
                    diatomz(zb+1)=KP_ratio*diatommix + (1-KP_ratio)*diatomz(zb+1);

                    POPmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*POPz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
                    POPz(1:zb)=POPmix;
                    POPz(zb+1)=KP_ratio*POPmix + (1-KP_ratio)*POPz(zb+1);

                    rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
                    TKE=0;
                    WmixIndicator=0;
                end %if (KP_ratio>=1)
            else
                WmixIndicator=0;
            end %if water column (not) already mixed
        end %while

        Aft_wind=sum(diff(rho)==0); %just a watch variable

    else % ice cover module
        XE_surf=(Tz(1)-Tf) * Cw * dz; %Daily heat accumulation into the first water layer (J m-2)
        Tz(1)=Tf; %Ensure that temperature of the first water layer is kept at freezing point
        TKE=0; %No energy for wind mixing under ice

        if (Wt(i,3)<Tf) %if air temperature is below freezing
            %Calculate ice surface temperature (Tice)
            if(WEQs==0) %if no snow
                alfa=1/(10*Hi);
                dHsi=0;
            else
                K_snow=2.22362*(rho_snow/1000)^1.885; %Yen (1981)
                alfa=(K_ice/K_snow)*(((rho_fw/rho_snow)*WEQs)/Hi);
                %Slush/snow ice formation (directly to ice)
                dHsi=max([0, Hi*(rho_ice/rho_fw-1)+WEQs]);
                Hsi=Hsi+dHsi;
            end
            Tice=(alfa*Tf+Wt(i,3))/(1+alfa);

            %Ice growth by Stefan's law
            Hi_new=sqrt((Hi+dHsi)^2+(2*K_ice/(rho_ice*L_ice))*(24*60*60)*(Tf-Tice));
            %snow fall
            dWEQnews=0.001*Wt(i,7); %mm->m
            dWEQs=dWEQnews-dHsi*(rho_ice/rho_fw); % new precipitation minus snow-to-snowice in snow water equivalent
            dHsi=0; %reset new snow ice formation
        else %if air temperature is NOT below freezing
            Tice=Tf;    %ice surface at freezing point
            dWEQnews=0; %No new snow
            if (WEQs>0)
                %snow melting in water equivalents
                dWEQs=-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_fw*L_ice)]);
                if ((WEQs+dWEQs)<0) %if more than all snow melts...
                    Hi_new=Hi+(WEQs+dWEQs)*(rho_fw/rho_ice); %...take the excess melting from ice thickness
                else
                    Hi_new=Hi; %ice does not melt until snow is melted away
                end
            else
                %total ice melting
                dWEQs=0;
                Hi_new=Hi-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_ice*L_ice)]);
                %snow ice part melting
                Hsi=Hsi-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_ice*L_ice)]);
                if (Hsi<=0)
                    Hsi=0;
                end
            end %if there is snow or not
        end %if air temperature is or isn't below freezing


        %Update ice and snow thicknesses
        Hi=Hi_new-(XE_surf/(rho_ice*L_ice)); %new ice thickness (minus melting due to heat flux from water)
        XE_surf=0; %reset energy flux from water to ice (J m-2 per day)
        WEQs=WEQs+dWEQs; %new snow water equivalent

        if(Hi<Hsi)
            Hsi=max(0,Hi);    %to ensure that snow ice thickness does not exceed ice thickness
            %(if e.g. much ice melting much from bottom)
        end


        if(WEQs<=0)
            WEQs=0; %excess melt energy already transferred to ice above
            rho_snow=rho_new_snow;
        else
            %Update snow density as weighed average of old and new snow densities
            rho_snow=rho_snow*(WEQs-dWEQnews)/WEQs + rho_new_snow*dWEQnews/WEQs;
            if (snow_compaction_switch==1)
                %snow compaction
                if (Wt(i,3)<Tf) %if air temperature is below freezing
                    rhos=1e-3*rho_snow; %from kg/m3 to g/cm3
                    delta_rhos=24*rhos*C1*(0.5*WEQs)*exp(-C2*rhos)*exp(-0.08*(Tf-0.5*(Tice+Wt(i,3))));
                    rho_snow=min([rho_snow+1e+3*delta_rhos, max_rho_snow]);  %from g/cm3 back to kg/m3
                else
                    rho_snow=max_rho_snow;
                end
            end
        end

        if(Hi<=0)
            IceIndicator=0;
            disp(['Ice-off, ' datestr(datenum(M_start)+i-1)])
            XE_melt=(-Hi-(WEQs*rho_fw/rho_ice))*rho_ice*L_ice/(24*60*60);
            %(W m-2) snow part is in case ice has melted from bottom leaving some snow on top (reducing XE_melt)
            Hi=0;
            WEQs=0;
            Tice=NaN;
            DoM(pp)=i;
            pp=pp+1;
        end

    end %of ice cover module

    %== P-partitioning in water==
    % TIPz=Pz + PPz; % Total inorg. phosphorus (excl. Chla and DOP) in the water column (mg m-3)
    % [Pz, trash]=Ppart(POCz./rho_sed,TIPz,Psat_L,Fmax_L,rho_sed,Fstable);
    % PPz=TIPz-Pz;

    %DIC-partitioning in water

    [CO2z,~] = carbonequilibrium(DICz,Tz,pH);

    % Relative dissolved oxygen concentration

    [O2_sat_rel, O2_sat_abs] = relative_oxygen(O2z,Tz,Wt(i,5),dz);

    %Initial freezing
    Supercooled=find(Tz<Tf);
    if (isempty(Supercooled)==0)
        %===NEW!!! (040707)
        if(Supercooled(1)~=1); disp('NOTE: non-surface subsurface supercooling'); end;
        InitIceEnergy=sum((Tf-Tz(Supercooled)).*Vz(Supercooled)*Cw);
        HFrazil=HFrazil+(InitIceEnergy/(rho_ice*L_ice))/Az(1);
        Tz(Supercooled)=Tf;

        if ((IceIndicator==0)&(HFrazil > Frazil2Ice_tresh))
            IceIndicator=1;
            Hi=Hi+HFrazil;
            HFrazil=0;
            DoF(qq)=i;
            disp(['Ice-on, ' datestr(datenum(M_start)+i-1)])
            qq=qq+1;
        end

        if (IceIndicator==1)
            Hi=Hi+HFrazil;
            HFrazil=0;
        end
        Tz(1)=Tf; %set temperature of the first layer to freezing point
        %======================

    end

    % Calculate pycnocline depth
    pycno_thres=0.1;  %treshold density gradient value (kg m-3 m-1)
    rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
    dRdz = [NaN; abs(diff(rho))];
    di=find((dRdz<(pycno_thres*dz)) | isnan(dRdz));
    %dRdz(di)=NaN;
    %TCz = nansum(zz .* dRdz) ./ nansum(dRdz);
    dRdz(di)=0; %modified for MATLAB version 7
    TCz = sum(zz .* dRdz) ./ sum(dRdz);

    %vector with S_res_epi above, and S_res_hypo below the pycnocline
    inx=find(zz <= TCz);
    S_resusp(inx)=S_res_epi;
    inx=find(zz > TCz);
    S_resusp(inx)=S_res_hypo;

    if (IceIndicator==1)
        S_resusp(:)=S_res_hypo;  %only hypolimnetic type of resuspension allowed under ice
    end

    if( isnan(TCz) & (IceIndicator==0) )
        S_resusp(:)=S_res_epi;   %if no pycnocline and open water, resuspension allowed from top to bottom
    end


    % WC chemistry:
    if any(isnan(O2z)) | any(isnan(Chlz)) | any(isnan(DOCz)) | any(isnan(NO3z)) | any(isnan(Fe3z)) | any(isnan(SO4z)) | any(isnan(NH4z)) | any(isnan(Fe2z)) | any(isnan(H2Sz)) | any(isnan(HSz)) | any(isnan(Pz)) | any(isnan(Al3z)) | any(isnan(PPz)) | any(isnan(Ca2z)) | any(isnan(CO2z))
        error('NaN')
    end

    if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
        error('NaN')
    end

    if wc_chemistry_module

        O2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(O2z, 31998.8);
        Chlz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Chlz, 30973.762);
        DOCz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(DOCz, 12010.7);
        NO3z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(NO3z, 62004);
        Fe3z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Fe3z, 106867.0); %Fe(OH)3
        SO4z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(SO4z, 96062);
        NH4z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(NH4z, 18038);
        Fe2z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Fe2z, 55845);
        H2Sz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(H2Sz, 34080.9);
        HSz     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(HSz, 33072.9);
        Pz      = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Pz, 30973.762); % PO4
        Al3z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Al3z, 78003.6); % Al(OH)3
        PPz     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(PPz, 30973.762); % PO4-s
        Ca2z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Ca2z, 80156.0); % Ca2+
        CO2z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(CO2z, 44009.5);
        DOPz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(DOPz, 30973.762); %DOPz
        Cz      = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Cz,  30973.762);
        POCz      = convert_mg_per_qubic_m_to_umol_per_qubic_cm(POCz,  12010.7);
        POPz      = convert_mg_per_qubic_m_to_umol_per_qubic_cm(POPz,  30973.762);

        % [Fe3z, Pz, PPz] = equilibrium_P_sorption(Fe3z, Pz, PPz, Kads);

        % Passing some MyLake results for chem module
        mylake_temp_results.Tz = Tz;
        mylake_temp_results.DayFrac = DayFrac;
        mylake_temp_results.lambdaz_wtot = lambdaz_wtot;
        mylake_temp_results.H_sw_z = H_sw_z;
        mylake_temp_results.H_sw_z_2 = H_sw_z_2;

        C0 = [O2z, Chlz, DOCz, NO3z, Fe3z, SO4z, NH4z, Fe2z, H2Sz, HSz, Pz, Al3z, PPz, Ca2z, CO2z, DOPz, Cz, POCz, POPz];

        [C_new, wc_rates_av] = wc_chemical_reactions_module(mylake_params, sediment_params, mylake_temp_results, C0, dt, sediment_params.n_of_time_steps_during_1_dt_of_myLake, wc_int_method);

        O2z  = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,1), 31998.8);
        Chlz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,2), 30973.762);
        DOCz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,3), 12010.7);
        NO3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,4), 62004);
        Fe3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,5), 106867.0); %Fe(OH)3
        SO4z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,6), 96062);
        NH4z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,7), 18038);
        Fe2z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,8), 55845);
        H2Sz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,9), 34080.9);
        HSz  = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,10), 33072.9);
        Pz   = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,11), 30973.762); % PO4
        Al3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,12), 78003.6);
        PPz  = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,13), 30973.762);
        Ca2z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,14), 80156.0); % Ca2+
        CO2z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,15), 44009.5);
        DOPz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,16), 30973.762); %DOPz
        Cz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,17), 30973.762);
        POCz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,18), 12010.7);
        POPz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,19), 30973.762);
    end

    if any(isnan(C_new))
        error('NaN')
    end

    if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
        error('NaN')
    end


    % sediment module
    if matsedlab_sediment_module
            % Making cells of params for using during coupling

            mylake_temp_results.Chlz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Chlz, 30973.762);
            mylake_temp_results.Cz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Cz,  30973.762);
            mylake_temp_results.POPz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(POPz, 30973.762);
            mylake_temp_results.DOCz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(DOCz, 12010.7);
            mylake_temp_results.DOPz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(DOPz, 30973.762);
            mylake_temp_results.O2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(O2z, 31998.8);
            mylake_temp_results.Pz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Pz, 30973.762);
            mylake_temp_results.Fe2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Fe2z, 55845);
            mylake_temp_results.NO3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(NO3z, 62004);
            mylake_temp_results.NH4z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(NH4z, 18038);
            mylake_temp_results.SO4z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(SO4z, 96062);
            mylake_temp_results.Fe3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Fe3z, 106867.0);
            mylake_temp_results.Ca2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Ca2z, 80156.0);
            mylake_temp_results.Al3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(Al3z, 78003.6);
            mylake_temp_results.PPz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(PPz, 30973.762);
            mylake_temp_results.POCz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(POCz,  12010.7);


        % Preparing units and estimate flux from [WC] ----> [Sediments]
        sediment_bc = update_sediment(mylake_temp_results, mylake_params, sediment_params);

        % Running sediment module
        [sediment_bioirrigation_fluxes, sediment_transport_fluxes, sediment_concentrations, sediment_additional_results] = sediment_v2(...
            sediment_concentrations, sediment_params, sediment_matrix_templates, sediment_bc);

        mylake_prev_results.Tz = Tz;
        mylake_prev_results.O2z = O2z;
        mylake_prev_results.Pz = Pz;
        mylake_prev_results.Fe2z = Fe2z;
        mylake_prev_results.NO3z = NO3z;
        mylake_prev_results.NH4z = NH4z;
        mylake_prev_results.DOPz = DOPz;
        mylake_prev_results.DOCz = DOCz;

        if any(isnan(O2z)) | any(isnan(Pz)) | any(isnan(Fe2z)) | any(isnan(NO3z)) | any(isnan(NH4z))
            error('NaN')
        end

        % Update WC:  [sediment] ----> [WC]
        [mylake_new_resutls] = update_wc(mylake_prev_results, sediment_concentrations, sediment_transport_fluxes, sediment_bioirrigation_fluxes, mylake_params, sediment_params);

        O2z = mylake_new_resutls.O2z;
        Pz = mylake_new_resutls.Pz;
        Fe2z = mylake_new_resutls.Fe2z;
        NO3z = mylake_new_resutls.NO3z;
        NH4z = mylake_new_resutls.NH4z;
        DOPz = mylake_new_resutls.DOPz;
        DOCz = mylake_new_resutls.DOCz;


        fields = fieldnames(sediment_concentrations);
        for fd_idx = 1:numel(fields)
            sediment_concentrations_zt.(fields{fd_idx})(:,i) = sediment_concentrations.(fields{fd_idx});
        end

        fields = fieldnames(sediment_transport_fluxes);
        for fd_idx = 1:numel(fields)
            sediment_transport_fluxes_zt.(fields{fd_idx})(:,i) = sediment_transport_fluxes.(fields{fd_idx});
        end

        fields = fieldnames(sediment_bioirrigation_fluxes);
        for fd_idx = 1:numel(fields)
            sediment_bioirrigation_fluxes_zt.(fields{fd_idx})(:,i) = sediment_bioirrigation_fluxes.(fields{fd_idx});
        end

        if mylake_params.rate_estimator_switch

            fields = fieldnames(sediment_additional_results.rates);
            for fd_idx = 1:numel(fields)
                sediment_additional_results_zt.rates.(fields{fd_idx})(:,i) = sediment_additional_results.rates.(fields{fd_idx});
            end

            fields = fieldnames(wc_rates_av);
            for fd_idx = 1:numel(fields)
                MyLake_results.rates.(fields{fd_idx})(:,i) = wc_rates_av.(fields{fd_idx});
            end
        else
            sediment_additional_results_zt.rates = false;
        end

    end

    % Output MyLake matrices
    Qst(:,i) = [Qsw Qlw Qsl]';
    Kzt(:,i) = [0;Kz];
    Tzt(:,i) = Tz;
    Czt(:,i) = Cz;
    POCzt(:,i) = POCz;
    Pzt(:,i) = Pz;
    Chlzt(:,i) = Chlz;
    PPzt(:,i) = PPz;
    DOPzt(:,i) = DOPz;
    DOCzt(:,i) = DOCz;
    DICzt(:,i) = DICz;
    O2zt(:,i) = O2z;

    NO3zt(:,i) = NO3z;
    NH4zt(:,i) = NH4z;
    SO4zt(:,i) = SO4z;
    HSzt(:,i) = HSz;
    H2Szt(:,i) = H2Sz;
    Fe2zt(:,i) = Fe2z;
    Ca2zt(:,i) = Ca2z;
    pHzt(:,i) = pHz;
    CH4zt(:,i) = CH4z;
    Fe3zt(:,i) = Fe3z;
    Al3zt(:,i) = Al3z;
    SiO4zt(:,i) = SiO4z;
    SiO2zt(:,i) = SiO2z;
    diatomzt(:,i) = diatomz;
    POPzt(:,i) = POPz;

    H_sw_zt(:,i) = diff([-H_sw_z; 0]);;
    H_sw_zt_2(:,i) = diff([-H_sw_z_2; 0]);;
    PAR_zt(:,i) = PAR_z;

    lvlDzt(:,i) = lvlD;

    % O2diffzt(:,i) = O2_diff;
    CO2zt(:,i) = CO2z;
    O2_sat_relt(:,i) = O2_sat_rel;
    O2_sat_abst(:,i) = O2_sat_abs;
    BODzt = 0; %for compatibility with the other code


    %Fokema
    %CDOMzt(:,i)=CDOMz;
    if (photobleaching==1)
        DOCzt1(:,i) = DOCz1_new; %Fokema-model DOC subpool 1
        DOCzt2(:,i) = DOCz2_new; %Fokema-model DOC subpool 2
        DOCzt3(:,i) = DOCz3_new; %Fokema-model DOC subpool 3
        DOC1tfrac(:,i) = DOC1frac; %Fokema-model subpool 1 fraction
        DOC2tfrac(:,i) = DOC2frac; %Fokema-model subpool 2 fraction
        DOC3tfrac(:,i) = DOC3frac; %Fokema-model subpool 3 fraction
        Daily_BB1t(:,i) = Daily_BB1; %Fokema-model subpool 1 daily bacterial decomposition
        Daily_BB2t(:,i) = Daily_BB2; %Fokema-model subpool 2 daily bacterial decomposition
        Daily_BB3t(:,i) = Daily_BB3; %Fokema-model subpool 3 daily bacterial decomposition
        Daily_PBt(:,i) = Daily_PB; %Fokema-model daily photobleaching
    end

    Qzt_sed(:,i) = Qz_sed./(60*60*24*dt); %(J m-2 day-1) -> (W m-2)
    lambdazt(:,i) = lambdaz_wtot_avg;

    surfaceflux(1,i) = surfflux; %Carbon dioxide surface flux
    % CO2_eqt(1,i) = CO2_eq;       %Carbon dioxide equilibrium concentration
    % K0t(:,i) = K0;               %Dissolved carbon doxide solubility coefficient
    % CO2_ppmt(:,i) = CO2_ppm;

    % O2fluxt(1,i) = O2flux;       %Oxygen surface flux
    % O2_eqt(1,i) = O2_eq;         %Oxygen equilibrium concentration
    % K0_O2t(1,i) = K0_O2;         %Dissolved oxygen solubility coefficient
    % dO2Chlt(:,i) = dO2_Chl;
    % dO2BODt(:,i) = dO2_BOD;
    % dO2SODt(:,i) = dO2_SOD;

    % testi1t(:,i) = O2_old;
    % testi2t(:,i) = O2_diff; testi3t(:,i) = O2_new;

    % P3zt_sed(:,i,1) = Pdz_store; %diss. P conc. in sediment pore water (mg m-3)
    % P3zt_sed(:,i,2) = Psz_store; %P conc. in inorganic sediment particles (mg kg-1 dry w.)
    % P3zt_sed(:,i,3) = Chlsz_store; %Chl conc. in organic sediment particles (mg kg-1 dry w.)
    % P3zt_sed(:,i,4) = F_IM; %VOLUME fraction of inorganic particles of total dry sediment
    % P3zt_sed(:,i,5) = Sedimentation_counter; %H_netsed_inorg; %Sedimentation (m/day) of inorganic particles of total dry sediment
    % P3zt_sed(:,i,6) = Resuspension_counter; %H_netsed_org; %Sedimentation (m/day) of organic particles of total dry sediment
    % P3zt_sed(:,i,7) = NewSedFrac; %(monitoring variables)

    % P3zt_sed_sc(:,i,1) = dPW_up; %(mg m-3 day-1) change in Pz due to exchange with pore water
    % P3zt_sed_sc(:,i,2) = dPP; %(mg m-3 day-1)
    % P3zt_sed_sc(:,i,3) = dChl_res; %(mg m-3 day-1)

    His(1,i) = Hi;
    His(2,i) = (rho_fw/rho_snow)*WEQs;
    His(3,i) = Hsi;
    His(4,i) = Tice;
    His(5,i) = Wt(i,3);
    His(6,i) = rho_snow;
    His(7,i) = IceIndicator;
    His(8,i) = HFrazil; %NEW!!!

    %Original MixStat matrix in v.1.2.1b

    % MixStat(1,i) = Iflw_S;
    % MixStat(2,i) = Iflw_TP;
    % MixStat(3,i) = sum(POCz.*Vz);
    %MixStat(4,i) = Growth_bioz(1);%mean(Growth_bioz(1:4)); %Obs! changed to apply to layers 1-4 only
    %MixStat(5,i) = Loss_bioz(1);%mean(Loss_bioz(1:4)); %Obs! changed to apply to layers 1-4 only
    %%MixStat(6,i) = Iflw;
    %MixStat(7:10,i) = NaN;

    % MixStat matrix from v.1.2 for figure output purposes

    MixStat(1,i) = Iflw_POC;
    MixStat(2,i) = Iflw_TP;
    MixStat(3,i) = lambdaz_wtot(2);%Iflw_DOC;
    Growth_bioz = 0; % NOTE: This part moved to reaction module
    Loss_bioz = 0; % NOTE: This part moved to reaction module
    MixStat(4,i) = mean(Growth_bioz); %Only for chlorophyll group 1 (a)
    MixStat(5,i) = mean(Loss_bioz);  %Only for chlorophyll group 1 (a)
    MixStat(6,i) = Iflw;
    if (IceIndicator == 1)
        MixStat(7:11,i) = NaN;
    else
        dum=interp1(zz,Pz,[0:0.1:4]);
        MixStat(7,i) = mean(dum); %diss-P conc. 0-4m in ice-free period

        dum=interp1(zz,Chlz,[0:0.1:4]);
        MixStat(8,i) = mean(dum); %Chla conc. 0-4m in ice-free period

        dum=interp1(zz,PPz,[0:0.1:4]);
        MixStat(9,i) = mean(dum); %particulate inorg. P conc. 0-4m in ice-free period

        dum=interp1(zz,DOPz,[0:0.1:4]);
        MixStat(10,i) = mean(dum); %dissolved organic P conc. 0-4m in ice-free period

        dum=interp1(zz,POCz,[0:0.1:4]);
        MixStat(11,i) = mean(dum); %particulate matter conc. 0-4m in ice-free period
    end

    MixStat(12,i) = TCz; %pycnocline depth

    MixStat(13,i) = 1e-6*Iflw*Iflw_TP; %total P inflow (kg day-1)
    %         if (Iflw>Vz(1))
    %             disp('Large inflow!!')
    %         end
    MixStat(14,i) = 1e-6*Iflw*(Pz(1)+PPz(1)+DOPz(1)+Chlz(1)+Cz(1)); %total P outflow (kg day-1)
    % MixStat(15,i) = sum(1e-6*Vz.*(delPP_inorg + delC_org)); %total P sink due to sedimentation (kg day-1)
    % MixStat(16,i) = sum(1e-6*(dPP+dPW_up).*Vz); %Internal P loading (kg day-1, excluding Chla)
    % MixStat(17,i) = sum(1e-6*dChl_res.*Vz); %Internal Chla loading (kg day-1)(resuspension 50/50 between the two groups)
    % MixStat(18,i)= sum(1e-6*Vz.*((Pz+PPz+DOPz+Chlz+Cz) - TP0)); %Net P change kg
    % MixStat(19,i)= sum(1e-6*((dPP+dPW_up-delPP_inorg+dChl_res-delC_org).*Vz - (1-F_sed_sld)*H_sed*(-diff([Az; 0])).*dPW_down)); %Net P flux from sediment kg
    % MixStat(20,i) = 1e-6*Iflw*(Iflw_TP-(Iflw_Chl+Iflw_C)./Y_cp-Iflw_DOP-Fstable*Iflw_S); %total algae-available P inflow (kg day-1)
    if (IceIndicator == 1)
        MixStat(21,i) = NaN;
    else
        dum=interp1(zz,Cz,[0:0.1:4]);
        MixStat(21,i) = mean(dum); %Chl group 2 conc. 0-4m in ice-free period
    end
    Growth_bioz_2 = 0;% NOTE: This part moved to reaction module
    Loss_bioz_2 = 0;  % NOTE: This part moved to reaction module
    MixStat(22,i) = mean(Growth_bioz_2); %For chlorophyll group 2
    MixStat(23,i) = mean(Loss_bioz_2);  %For chlorophyll group 2

end; %for i = 1:length(tt)

stamp = datetime('now');

%Saving sediment values
if matsedlab_sediment_module           % MATSEDLAB sediment module
    sediment_results.concentrations = sediment_concentrations_zt;
    sediment_results.z = sediment_params.x';
    sediment_results.Bioirrigation_fx_zt = sediment_bioirrigation_fluxes_zt;
    sediment_results.mylake_params = mylake_params;
    sediment_results.params = sediment_params;
    sediment_results.sediment_transport_fluxes = sediment_transport_fluxes_zt;
    sediment_results.days = datenum(M_start):datenum(M_stop);
    sediment_results.m_start = M_start;
    sediment_results.m_stop = M_stop;
    sediment_results.rates = sediment_additional_results_zt.rates;
    sediment_results.date_of_run = stamp;
else

    sediment_results = {};

end

d_O2zt = (diff(O2zt'))';
int_R_O2dz = integrate_over_depth(d_O2zt, dz);


MyLake_results.Qst = Qst;
MyLake_results.K = Kzt;
MyLake_results.T = Tzt;
MyLake_results.concentrations.P = Pzt;
MyLake_results.concentrations.PP = PPzt;
MyLake_results.concentrations.C = Czt;
MyLake_results.concentrations.Chl = Chlzt;
MyLake_results.concentrations.POP = POPzt;
MyLake_results.concentrations.DOP = DOPzt;
MyLake_results.concentrations.DOC = DOCzt;
MyLake_results.concentrations.DIC = DICzt;
MyLake_results.concentrations.CO2 = CO2zt;
MyLake_results.concentrations.O2 = O2zt;
MyLake_results.concentrations.NO3 = NO3zt;
MyLake_results.concentrations.NH4 = NH4zt;
MyLake_results.concentrations.Fe3 = Fe3zt;
MyLake_results.concentrations.Fe2 = Fe2zt;
MyLake_results.concentrations.SO4 = SO4zt;
MyLake_results.concentrations.HS = HSzt;
MyLake_results.concentrations.H2S = H2Szt;
MyLake_results.concentrations.CH4 = CH4zt;
MyLake_results.concentrations.Ca2 = Ca2zt;
MyLake_results.concentrations.pH = pHzt;
MyLake_results.concentrations.POC = POCzt;
MyLake_results.concentrations.Al3 = Al3zt;
MyLake_results.concentrations.SiO4 = SiO4zt;
MyLake_results.concentrations.SiO2 = SiO2zt;
MyLake_results.concentrations.diatom = diatomzt;
MyLake_results.d_O2zt = d_O2zt;
MyLake_results.int_R_O2dz = int_R_O2dz;
MyLake_results.O2_diffzt = O2_diffzt;
MyLake_results.O2_sat_relt = O2_sat_relt;
MyLake_results.O2_sat_abst = O2_sat_abst;
MyLake_results.Qzt_sed = Qzt_sed;
MyLake_results.lambdazt = lambdazt;
MyLake_results.P3zt_sed = P3zt_sed;
MyLake_results.P3zt_sed_sc = P3zt_sed_sc;
MyLake_results.His = His;
MyLake_results.MixStat = MixStat;
MyLake_results.H_sw = H_sw_zt;
MyLake_results.H_sw_2 = H_sw_zt_2;
MyLake_results.PAR_zt = PAR_zt;
MyLake_results.lvlDzt = lvlDzt;
MyLake_results.CDOMzt = CDOMzt;
MyLake_results.DOCzt1 = DOCzt1;
MyLake_results.DOCzt2 = DOCzt2;
MyLake_results.DOCzt3 = DOCzt3;
MyLake_results.DOC1tfrac = DOC1tfrac;
MyLake_results.DOC2tfrac = DOC2tfrac;
MyLake_results.DOC3tfrac = DOC3tfrac;
MyLake_results.Daily_BB1t = Daily_BB1t;
MyLake_results.Daily_BB2t = Daily_BB2t;
MyLake_results.Daily_BB3t = Daily_BB3t;
MyLake_results.Daily_PBt = Daily_PBt;
MyLake_results.surfaceflux = surfaceflux;
MyLake_results.CO2_eqt = CO2_eqt;
MyLake_results.CO2_ppmt = CO2_ppmt;
MyLake_results.K0t = K0t;
MyLake_results.O2fluxt = O2fluxt;
MyLake_results.O2_eqt = O2_eqt;
MyLake_results.K0_O2t = K0_O2t;
MyLake_results.dO2Chlt = dO2Chlt;
MyLake_results.dO2BODt = dO2BODt;
MyLake_results.dfloc_DOC = dfloc_DOC;
MyLake_results.testi1t = testi1t;
MyLake_results.testi2t = testi2t;
MyLake_results.z = zz;
MyLake_results.days = datenum(M_start):datenum(M_stop);
MyLake_results.params = mylake_params;
MyLake_results.params.Phys_par = Phys_par;
MyLake_results.params.Bio_par = Bio_par;
MyLake_results.Inflw = Inflw;
MyLake_results.Wt = Wt;
MyLake_results.m_start = M_start;
MyLake_results.m_stop = M_stop;
MyLake_results.date_of_run = stamp;

runtime=toc;

%disp(['Total model runtime: ' int2str(floor(runtime/60)) ' min ' int2str(round(mod(runtime,60))) ' s']);
%disp(['Reduced SS load due to inconsistencies: '  num2str(round(SS_decr)) ' kg']);

% >>>>>> End of the time loop >>>>>>

% Below are the two functions for calculating tridiagonal matrix Fi for solving the
% 1) diffusion equation (tridiag_DIF_v11), and
% 2) advection-diffusion equation (tridiag_HAD_v11) by fully implicit hybrid exponential numerical scheme,
% based on Dhamotharan et al. 1981,
%'Unsteady one-dimensional settling of suspended sediments', Water Resources Research 17(4), 1125-1132
% code checked by TSA, 16.03.2004


%Inputs:
% Kz    diffusion coefficient at layer interfaces (plus surface) N (N,1)
% U     vertical settling velocity (scalar)
% Vz    layer volumes (N,1)
% Az    layer interface areas (N,1)
% dz    grid size
% dt    time step

%Output:
% Fi    tridiagonal matrix for solving new profile Cz

% az = (dt/dz) * [0; Kz] .* (Az ./ Vz);
% bz = (dt/dz) * [Kz; 0] .* ([Az(2:end); 0] ./ Vz);
% Gi = [-bz (1 + az + bz) -az];

%=== DIFFUSIVE EQUATION ===
function [Fi] = tridiag_DIF_v11(Kz,Vz,Az,dz,dt)

Nz=length(Vz); %number of grid points/layers

% Linearized heat conservation equation matrix (diffusion only)
az = (dt/dz) * Kz .* (Az ./ Vz);                                        %coefficient for i-1
cz = (dt/dz) * [Kz(2:end); NaN] .* ([Az(2:end); NaN] ./ Vz);            %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i+1
%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged
bz(1)= 1 + az(1) + cz(1);


%Boundary conditions, bottom

%az(end) remains unchanged
cz(end) = 0;
bz(end) = 1 + az(end) + cz(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function


%=== ADVECTIVE-DIFFUSIVE EQUATION ===
function [Fi] = tridiag_HAD_v11(Kz,U,Vz,Az,dz,dt)

if (U<0)
    error('only positive (downward) velocities allowed')
end

if (U==0)
    U=eps; %set Vz next to nothing (=2.2204e-016) in order to avoid division by zero
end

Nz=length(Vz); %number of grid points/layers

theta=U*(dt/dz);

az = theta.*(1 + (1./(exp( (U*Vz)./(Kz.*Az) ) - 1)));                   %coefficient for i-1
cz = theta./(exp( (U*Vz)./([Kz(2:end); NaN].*[Az(2:end); NaN]) ) - 1);  %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i

%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged
bz(1) = 1 + theta + cz(1);

%Boundary conditions, bottom

%az(end) remains unchanged
cz(end) = 0;
bz(end) = 1 + az(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function



function [C_new, rates] = wc_chemical_reactions_module(mylake_params, sediment_params, mylake_temp_results, C0, dt_mylake, n_of_time_steps_during_1_dt_of_myLake, method)
    % ts - how many time steps during 1 day
    % dt - time step in chemical module [years]
    % dt_mylake - time step in MyLake [days]
    % n_of_time_steps_during_1_dt_of_myLake - amount of steps during of 1 dt of MyLake;


    dt = (dt_mylake/365) / n_of_time_steps_during_1_dt_of_myLake;

    if method == 0
        [C_new, rates] = rk4(mylake_params, sediment_params, mylake_temp_results, C0, dt, n_of_time_steps_during_1_dt_of_myLake);
    elseif method == 1
        [C_new, rates] = butcher5(mylake_params, sediment_params, mylake_temp_results, C0, dt, n_of_time_steps_during_1_dt_of_myLake);
    end
    C_new = (C_new>0).*C_new;

%end of function


%% rk4: Runge-Kutta 4th order integration
function [C_new, rates_av] = rk4(mylake_params, sediment_params, mylake_temp_results, C0, dt,n)
    % ts - time step in chemical module [years]
    % dt - time step in MyLake [days]
    % (1/365/ts) = is how many steps during 1 dt of Mylake
    % (dt/365) = is conversion of [days] to [years]

    for i = 1:n
        [dcdt_1, r_1] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0, dt);
        k_1 = dt.*dcdt_1;
        [dcdt_2, r_2] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0+0.5.*k_1, dt);
        k_2 = dt.*dcdt_2;
        [dcdt_3, r_3] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0+0.5.*k_2, dt);
        k_3 = dt.*dcdt_3;
        [dcdt_4, r_4] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0+k_3, dt);
        k_4 = dt.*dcdt_4;
        C_new = C0 + (k_1+2.*k_2+2.*k_3+k_4)/6;
        C0 = C_new;

        % average rate
        if mylake_params.rate_estimator_switch
            fields = fieldnames(r_1);
            for fld_idx = 1:numel(fields)
              r.(fields{fld_idx}) = (r_1.(fields{fld_idx}) + 2*r_2.(fields{fld_idx}) + 2*r_3.(fields{fld_idx}) + r_4.(fields{fld_idx}))/6;
            end

            rates(i) = r;
        end
    end

    if mylake_params.rate_estimator_switch
        fields = fieldnames(rates);
        for i = 1:numel(fields)
          rates_av.(fields{i}) = 0;
          for j=1:ts-1
              rates_av.(fields{i}) = rates_av.(fields{i}) + rates(j).(fields{i});
          end
          rates_av.(fields{i}) = rates_av.(fields{i})/(ts-1);
        end
    else
        rates_av = false;
    end

%% butcher5: Butcher's Fifth-Order Runge-Kutta
function [C_new, rates_av] = butcher5(mylake_params, sediment_params, mylake_temp_results, C0, dt,n)

    for i = 1:n
        [dcdt_1, r_1] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0, dt);
        k_1 = dt.*dcdt_1;
        [dcdt_2, r_2] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 + 1/4.*k_1, dt);
        k_2 = dt.*dcdt_2;
        [dcdt_3, r_3] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 + 1/8.*k_1 + 1/8.*k_2, dt);
        k_3 = dt.*dcdt_3;
        [dcdt_4, r_4] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 - 1/2.*k_2 + k_3, dt);
        k_4 = dt.*dcdt_4;
        [dcdt_5, r_5] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 + 3/16.*k_1 + 9/16.*k_4, dt);
        k_5 = dt.*dcdt_5;
        [dcdt_6, r_6] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 - 3/7.*k_1 + 2/7.*k_2 + 12/7.*k_3 - 12/7.*k_4 + 8/7.*k_5, dt);
        k_6 = dt.*dcdt_6;
        C_new = C0 + (7.*k_1 + 32.*k_3 + 12.*k_4 + 32.*k_5 + 7.*k_6)/90;
        C0 = C_new;

        if mylake_params.rate_estimator_switch
            fields = fieldnames(r_1);
            for fld_idx = 1:numel(fields)
              r.(fields{fld_idx}) = (7*r_1.(fields{fld_idx}) + 32*r_3.(fields{fld_idx}) + 12*r_4.(fields{fld_idx}) + 32*r_5.(fields{fld_idx}) + 7*r_6.(fields{fld_idx}))/90;
            end
            rates(i) = r;
        end
    end

    if mylake_params.rate_estimator_switch
        fields = fieldnames(rates);
        for i = 1:numel(fields)
          rates_av.(fields{i}) = 0;
          for j=1:ts-1
              rates_av.(fields{i}) = rates_av.(fields{i}) + rates(j).(fields{i});
          end
          rates_av.(fields{i}) = rates_av.(fields{i})/(ts-1);
        end
    else
        rates_av = false;
    end


function [int_rate] = integrate_over_depth(R, dz)
%% integrate_over_depth: integrates the rates of reaction over the depth
% R - rate of interest
% dz - the mesh size
int_rate = sum(R)*dz;
%end of function

function C = convert_mg_per_qubic_m_to_umol_per_qubic_cm(C,M_C)
    C = C./M_C;
%end of function

function C = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C,M_C)
    C = C.*M_C;
%end of function

function [Fe3z, Pz, PPz] = equilibrium_P_sorption(Fe3z, Pz, PPz, Kads)
    x0 = [PPz(1), Pz(1)];
    for i=1:length(Fe3z)
        Stot(i) = Fe3z(i)+PPz(i);
        Atot(i) = Pz(i) + PPz(i);
        f = @(x)p_langmuir(x, Stot(i), Atot(i), Kads);

        options = optimoptions('fsolve','Display','none', 'MaxIterations', 1000);
        x = fsolve(f,x0, options);
        PPz(i) = x(1);
        Pz(i) = x(2);
        Fe3(i) = Stot(i) - PPz(i);
        x0=x;
    end

function F = p_langmuir(x, Stot, Atot, Kads)
    SA = x(1);
    A = x(2);
    F(1) = Stot.*Kads.*A/(1+Kads.*A) - SA;
    F(2) = A + SA - Atot;

function F = P_equlibrium_logK(x, param, Sum_Fe, Sum_P)
    pH = param;
    Fe3_H2P = x(1);
    Fe3_HP = x(2);
    Fe3_P = x(3);
    Fe3 = x(4);
    P = x(5);
    F(1) = log10(Fe3_H2P) - log10(Fe3) - log10(P) - 31.29 + 3 * pH;
    F(2) = log10(Fe3_HP) - log10(Fe3) - log10(P) - 25.39 + 2 * pH;
    F(3) = log10(Fe3_P) - log10(Fe3) - log10(P) - 17.72 +  pH  ;
    F(4) = Fe3_P + Fe3_HP + Fe3_H2P + Fe3 - Sum_Fe;
    F(5) = Fe3_P + Fe3_HP + Fe3_H2P + P - Sum_P;

function [dcdt, r] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C, dt)
% parameters for water-column chemistry
% NOTE: the rates are the same as in sediments (per year!! not day) except microbial which are factorized due to lower concertation of bacteria in WC then in sediments. Units are per "year" due to time step is in year units too;

    if any(isnan(C))
        error('NaN')
    end

    dcdt=zeros(size(C));

    O2z = C(:,1) .* (C(:,1)>0);
    Chlz = C(:,2) .* (C(:,2)>0);
    DOCz = C(:,3) .* (C(:,3)>0); % NOTE: POP: DOC 1:20 ratio. Conversion happens at SWI flux
    NO3z = C(:,4) .* (C(:,4)>0);
    Fe3z = C(:,5) .* (C(:,5)>0);
    SO4z = C(:,6) .* (C(:,6)>0);
    NH4z = C(:,7) .* (C(:,7)>0);
    Fe2z = C(:,8) .* (C(:,8)>0);
    H2Sz = C(:,9) .* (C(:,9)>0);
    HSz = C(:,10) .* (C(:,10)>0);
    Pz = C(:,11) .* (C(:,11)>0);
    Al3z = C(:,12) .* (C(:,12)>0);
    PPz = C(:,13) .* (C(:,13)>0);
    Ca2z = C(:,14) .* (C(:,14)>0);
    CO2z = C(:,15) .* (C(:,15)>0);
    DOPz = C(:,16) .* (C(:,16)>0);
    Cz = C(:,17) .* (C(:,17)>0);
    POCz = C(:,18) .* (C(:,18)>0);
    POPz = C(:,19) .* (C(:,19)>0);

    % ============ water-column module ============
    Km_O2         = sediment_params.Km_O2;
    Km_NO3        = sediment_params.Km_NO3;
    Km_FeOH3      = sediment_params.Km_FeOH3;
    Km_FeOOH      = sediment_params.Km_FeOOH;
    Km_SO4        = sediment_params.Km_SO4;
    Km_oxao       = sediment_params.Km_oxao;
    Km_amao       = sediment_params.Km_amao;
    Kin_O2        = sediment_params.Kin_O2;
    Kin_NO3       = sediment_params.Kin_NO3;
    Kin_FeOH3     = sediment_params.Kin_FeOH3;
    k_amox        = sediment_params.k_amox;
    k_Feox        = sediment_params.k_Feox;
    k_Sdis        = sediment_params.k_Sdis;
    k_Spre        = sediment_params.k_Spre;
    k_alum        = sediment_params.k_alum;
    k_pdesorb_c   = sediment_params.k_pdesorb_c;
    k_pdesorb_a   = sediment_params.k_pdesorb_a;
    k_pdesorb_b   = sediment_params.k_pdesorb_b;
    k_rhom        = sediment_params.k_rhom;
    k_tS_Fe       = sediment_params.k_tS_Fe;
    Ks_FeS        = sediment_params.Ks_FeS;
    k_Fe_dis      = sediment_params.k_Fe_dis;
    k_Fe_pre      = sediment_params.k_Fe_pre;
    k_apa         = sediment_params.k_apa;
    kapa          = sediment_params.kapa;
    k_oms         = sediment_params.k_oms;
    k_tsox        = sediment_params.k_tsox;
    k_FeSpre      = sediment_params.k_FeSpre;
    accel         = sediment_params.accel;
    f_pfe         = sediment_params.f_pfe;
    Cx1           = sediment_params.Cx1;
    Ny1           = sediment_params.Ny1;
    Pz1           = sediment_params.Pz1;
    Cx2           = sediment_params.Cx2;
    Ny2           = sediment_params.Ny2;
    Pz2           = sediment_params.Pz2;
    Cx3           = sediment_params.Cx3;
    Ny3           = sediment_params.Ny3;
    Pz3           = sediment_params.Pz3;


    T_ref = mylake_params.T_ref;
    Q10 = mylake_params.Q10;
    dop_twty = mylake_params.dop_twty;
    g_twty = mylake_params.g_twty;
    m_twty = mylake_params.m_twty;
    g_twty_2 = mylake_params.g_twty_2;
    m_twty_2 = mylake_params.m_twty_2;
    P_half = mylake_params.P_half;
    P_half_2 = mylake_params.P_half_2;
    theta_m = mylake_params.theta_m;
    dz = mylake_params.dz;
    floculation_switch = mylake_params.floculation_switch;

    Tz = mylake_temp_results.Tz;
    DayFrac = mylake_temp_results.DayFrac;
    lambdaz_wtot = mylake_temp_results.lambdaz_wtot;
    H_sw_z = mylake_temp_results.H_sw_z;
    H_sw_z_2 = mylake_temp_results.H_sw_z_2;


    k_POP         = sediment_params.k_POP * mylake_params.wc_factor;
    k_POC         = sediment_params.k_POC * mylake_params.wc_factor;
    k_DOP        = sediment_params.k_DOP * mylake_params.wc_factor;
    k_DOC        = sediment_params.k_DOC * mylake_params.wc_factor;
    k_POP_q10       = k_POP .* Q10.^((Tz-T_ref)/10);
    k_POC_q10       = k_POC .* Q10.^((Tz-T_ref)/10);
    k_DOP_q10       = k_DOP .* Q10.^((Tz-T_ref)/10);
    k_DOC_q10       = k_DOC .* Q10.^((Tz-T_ref)/10);


    % MyLake "old" chemistry:
    % Conversion of units to "per year" and umoles
    dop_twty_y = dop_twty*365; % In MyLake parameters with set to 0;
    g_twty_y = g_twty*365;
    m_twty_y = m_twty*365;
    g_twty_y_2 = g_twty_2*365;
    m_twty_y_2 = m_twty_2*365;
    P_half_molar = P_half/30973.762;
    P_half_2_molar = P_half_2/30973.762;

    % DOP: ( we treed this as bio reaction, see rate "Rc")
    R_dDOP = 0; % dop_twty_y .* DOPz .* theta_m.^(Tz-20);  %Mineralisation to P

    % Chlz:
    Growth_bioz=g_twty_y*theta_m.^(Tz-20) .* (Pz./(P_half_molar+Pz)) .* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z; 0]);
    Loss_bioz=m_twty_y*theta_m.^(Tz-20);
    R_bioz = Growth_bioz-Loss_bioz;


    % Cz:
    Growth_bioz_2=g_twty_y_2*theta_m.^(Tz-20) .* (Pz./(P_half_2_molar+Pz)) .* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z_2; 0]);
    Loss_bioz_2=m_twty_y_2*theta_m.^(Tz-20);
    R_bioz_2 = Growth_bioz_2-Loss_bioz_2;

    %flocculation
    if (floculation_switch==1) %Fokema
        dfloc_DOC = 0.030 .* DOCz;
        dfloc_DOP = 0.030 .* DOPz;
    else
        dfloc_DOC = 0;
        dfloc_DOP = 0;
    end


    R_dChl_growth =  Chlz .* R_bioz; %Chl a growth source
    R_dCz_growth =  Cz .* R_bioz_2;

    %Oxygen production in phytoplankton growth
    R_dO2_Chl = R_dChl_growth+R_dCz_growth;


    % New chemistry

    tot_FeOH3 = PPz + Fe3z;

    f_O2    = O2z ./  (Km_O2 + O2z) ;
    f_NO3   = NO3z ./  (Km_NO3 + NO3z) .* Kin_O2 ./ (Kin_O2 + O2z) ;
    f_FeOH3 = tot_FeOH3 ./  (Km_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3z) .* Kin_O2 ./ (Kin_O2 + O2z) ;
    f_SO4 = SO4z ./ (Km_SO4 + SO4z) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3z) .* Kin_O2 ./ (Kin_O2 + O2z);

    Sum_H2S = H2Sz + HSz;

    part_PO4ads_tot_Fe = PPz ./ (tot_FeOH3+1e-16); % Avoid division by zero

    R1a =  0; % k_POP_q10  .* Chlz .* f_O2 * accel;
    R1b =  0; % k_POP_q10  .* Cz .* f_O2 * accel;
    R1c =  k_DOP_q10  .* DOPz .* f_O2 .* accel;
    R1d =  k_DOC_q10 .* DOCz .* f_O2 .* accel;
    R1e =  k_POP_q10 .* POPz .* f_O2 .* accel;
    R1f =  k_POC_q10 .* POCz .* f_O2 .* accel;

    R2a =  0; % k_POP_q10  .* Chlz .* f_NO3 .* accel;
    R2b =  0; % k_POP_q10  .* Cz .* f_NO3 .* accel;
    R2c =  k_DOP_q10  .* DOPz .* f_NO3 .* accel;
    R2d =  k_DOC_q10 .* DOCz .* f_NO3 .* accel;
    R2e =  k_POP_q10 .* POPz .* f_NO3 .* accel;
    R2f =  k_POC_q10 .* POCz .* f_NO3 .* accel;

    R3a_Fe =  (1 - part_PO4ads_tot_Fe) .* 0; % k_POP_q10  .* Chlz .* f_FeOH3;
    R3b_Fe =  (1 - part_PO4ads_tot_Fe) .* 0; % k_POP_q10  .* Cz .* f_FeOH3;
    R3c_Fe =  (1 - part_PO4ads_tot_Fe) .* k_DOP_q10  .* DOPz .* f_FeOH3;
    R3d_Fe =  (1 - part_PO4ads_tot_Fe) .* k_DOC_q10 .* DOCz .* f_FeOH3;
    R3e_Fe =  (1 - part_PO4ads_tot_Fe) .* k_POP_q10 .* POPz .* f_FeOH3;
    R3f_Fe =  (1 - part_PO4ads_tot_Fe) .* k_POC_q10 .* POCz .* f_FeOH3;
    R3a_P =  part_PO4ads_tot_Fe .* 0; % k_POP_q10  .* Chlz .* f_FeOH3;
    R3b_P =  part_PO4ads_tot_Fe .* 0; % k_POP_q10  .* Cz .* f_FeOH3;
    R3c_P =  part_PO4ads_tot_Fe .* k_DOP_q10  .* DOPz .* f_FeOH3;
    R3d_P =  part_PO4ads_tot_Fe .* k_DOC_q10 .* DOCz .* f_FeOH3;
    R3e_P =  part_PO4ads_tot_Fe .* k_POP_q10 .* POPz .* f_FeOH3;
    R3f_P =  part_PO4ads_tot_Fe .* k_POC_q10 .* POCz .* f_FeOH3;
    R3a = R3a_P + R3a_Fe;
    R3b = R3b_P + R3b_Fe;
    R3c = R3c_P + R3c_Fe;
    R3d = R3d_P + R3d_Fe;
    R3e = R3e_P + R3e_Fe;
    R3f = R3f_P + R3f_Fe;

    R5a =  0; % k_POP_q10  .* Chlz .* f_SO4;
    R5b =  0; % k_POP_q10  .* Cz .* f_SO4;
    R5c =  k_DOP_q10  .* DOPz .* f_SO4;
    R5d =  k_DOC_q10 .* DOCz .* f_SO4;
    R5e =  k_POP_q10 .* POPz .* f_SO4;
    R5f =  k_POC_q10 .* POCz .* f_SO4;

    Ra  = R1a+R2a+R3a+R5a;
    Rb  = R1b+R2b+R3b+R5b;
    Rc  = R1c+R2c+R3c+R5c;
    Rd  = R1d+R2d+R3d+R5d;
    Re  = R1e+R2e+R3e+R5e;
    Rf  = R1f+R2f+R3f+R5f;

    R1  = R1a+R1b+R1c+R1d+R1e+R1f;
    R2  = R2a+R2b+R2c+R2d+R2e+R2f;
    R3  = R3a+R3b+R3c+R3d+R3e+R3f;
    R5  = R5a+R5b+R5c+R5d+R5e+R5f;
    R11  = k_tsox .* O2z .* Sum_H2S;
    R12  = k_tS_Fe .* Fe3z .* Sum_H2S;
    R13  = k_Feox .* Fe2z .* O2z;
    % NOTE: Due to the reaction is too fast and could cause overshooting:
    % we need to make this check if R*dt > Conc of source:
    % R13 = (R13.*dt < Fe2z).*R13 + (R13.*dt > Fe2z).* Fe2z ./ (dt) * 0.5;
    % R13 = (R13.*dt < O2z).*R13 + (R13.*dt > O2z).* O2z ./ (dt) * 0.5;

    % R14  = k_amox .* O2z ./ (Km_oxao + O2z) .* NH4z ./ (Km_amao + NH4z);   % NOTE: Doesnt work - Highly unstable.
    R14 = k_amox  .* O2z .* NH4z;
    % R14 = (R14.*dt < NH4z).*R14 + (R14.*dt > NH4z).* NH4z ./ (dt) * 0.5;
    % R14 = (R14.*dt < O2z).*R14 + (R14.*dt > O2z).* O2z ./ (dt) * 0.5;

    R21a = 0; %k_oms .* Sum_H2S .* Chlz;
    R21b = 0; %k_oms .* Sum_H2S .* Cz;
    R21c = k_oms .* Sum_H2S .* DOPz;
    R21d = k_oms .* Sum_H2S .* DOCz;
    R21e = k_oms .* Sum_H2S .* POPz;
    R21f = k_oms .* Sum_H2S .* POCz;

    R23  = 0; % NOTE: no FeS
    R24  = 0; % NOTE: no FeS
    R25a = 0; % NOTE: no FeS
    R25b = 0; % NOTE: no FeS

    R31a = k_pdesorb_a .* Fe3z .* Pz;
    % R31b = f_pfe .* (4 * R3 + 2 * R12);
    R31b = 4*(Cx1*R3a_P + Cx1*R3b_P + Cx1*R3c_P + Cx2*R3d_P+ Cx2*R3e_P+ Cx3*R3f_P);

    % R31b = (R31b.*dt < PPz).*R31b + (R31b.*dt > PPz).* PPz ./ (dt) * 0.5;

    R32a = 0; % No FeOOH in WC
    R32b = 0; % No FeOOH in WC
    R33a = k_pdesorb_c .* Pz .* Al3z; % NOTE: No separate pool for sorbed P on aluminum in WC
    R33a = 0; % NOTE: No separate pool for sorbed P on aluminum in WC
    R33b = 0; % NOTE: the rate is unknown
    R34  = k_apa .* (Pz - kapa); % NOTE: no Ca3PO42 pool in WC
    R34  = (R34 >= 0) .* R34;




    % saving rates
    r.R1a = R1a; r.R1b = R1b; r.R1c = R1c; r.R1d = R1d; r.R1e = R1e; r.R1f = R1f; r.R2a = R2a; r.R2b = R2b; r.R2c = R2c; r.R2d = R2d; r.R2e = R2e; r.R2f = R2f; r.R3a = R3a; r.R3b = R3b; r.R3c = R3c; r.R3d = R3d; r.R3e = R3e; r.R3f = R3f; r.R5a = R5a; r.R5b = R5b; r.R5c = R5c; r.R5d = R5d; r.R5e = R5e; r.R5f = R5f; r.Ra = Ra; r.Rb = Rb; r.Rc = Rc; r.Rd = Rd; r.Re = Re; r.Rf = Rf; r.R1 = R1; r.R2 = R2; r.R3 = R3; r.R5 = R5; r.R11 = R11; r.R12 = R12; r.R13  = R13; r.R14 = R14; r.R21a = R21a; r.R21b = R21b; r.R21c = R21c; r.R21d = R21d;; r.R21e = R21e;; r.R21f = R21f; r.R23 = R23; r.R24 = R24; r.R25a = R25a; r.R25b  = R25b; r.R31a = R31a; r.R31b  = R31b; r.R32a = R32a; r.R32b = R32b; r.R33a = R33a; r.R33b = R33b; r.R34 = R34;


    dcdt(:,1)  = -0.25*R13  - R11 - 2*R14 - (Cx1*R1a + Cx1*R1b + Cx2*R1c + Cx3*R1d+ Cx2*R1e+ Cx3*R1f) - 3*R23 + Cx1 * R_dO2_Chl; % O2z
    dcdt(:,2)  = -Ra - R21a + R_dChl_growth;% Chlz
    dcdt(:,3)  = -Rd - R21d - dfloc_DOC;% DOCz
    dcdt(:,4)  = - 0.8*(Cx1*R2a + Cx1*R2b + Cx2*R2c + Cx3*R2d+ Cx2*R2e + Cx3*R2f) + R14 - Ny1 * (R_dChl_growth + R_dCz_growth); % NO3z
    dcdt(:,5)  = - 4*(Cx1*R3a_Fe + Cx1*R3b_Fe + Cx2*R3c_Fe + Cx3*R3d_Fe+ Cx2*R3e_Fe+ Cx3*R3f_Fe) - 2*R12  + R13 - R31a; % Fe3z
    dcdt(:,6)  = - 0.5*(Cx1*R5a + Cx1*R5b + Cx2*R5c + Cx3*R5d+ Cx2*R5e+ Cx3*R5f) + R11 ; % SO4z
    dcdt(:,7)  =  (Ny1 * Ra + Ny1 * Rb + Ny2 * Rc + Ny3 * Rd+ Ny2 * Re+ Ny3 * Rf) - R14 ;% NH4z
    dcdt(:,8)  = 4*(Cx1*R3a + Cx1*R3b + Cx2*R3c + Cx3*R3d + Cx2*R3e+ Cx3*R3f) + 2*R12 - R13 + R25b - R25a; % Fe2z
    dcdt(:,9)  =  0;% H2Sz
    dcdt(:,10) = 0.5*(Cx1*R5a + Cx1*R5b + Cx2*R5c + Cx3*R5d+ Cx2*R5e+ Cx3*R5f) - R11 - R12  - R21a - R21b - R21c - R21d - R21e - R21f + R25b - R25a - R24 ;% HSz
    dcdt(:,11) = (Pz1 * Ra + Pz1 * Rb + Pz2 * Rc + Pz3 * Rd+ Pz2 * Re+ Pz3 * Rf) - R33a + R33b - R31a - R32a + R31b + R32b - 2*R34 + R_dDOP - Pz1 * (R_dChl_growth + R_dCz_growth);% Pz
    dcdt(:,12) = -R33a ;% Al3z
    dcdt(:,13) = R31a - R31b;% PPz
    dcdt(:,14) = -3*R34 ;% Ca2z
    dcdt(:,15) = -Rf - R21f;% CO2z
    dcdt(:,16) = -Rc - R21c - R_dDOP - dfloc_DOP;% DOPz
    dcdt(:,17) = -Rb - R21b + R_dCz_growth;% Cz
    dcdt(:,18) = dfloc_DOC - Rf; % POCz
    dcdt(:,19) = - Re + dfloc_DOP;% POPz
%end of function
