function[O2_sat_rel, O2_sat_abs] = relative_oxygen(O2z,Tz,Pa,dz)
% Calculates the relative and absolute saturations of dissolved oxygen

%Inputs
%O2z         O2 concentration (mg/m^3)
%Pa          Air pressure (mbar = hPa)
%Tz          Water temperature (oC)
%dz          Grid step (m)

%Outputs
%O2_sat_rel  Relative oxygen saturation
%O2_sat_abs  Absolute oxygen saturation


global ies80;
%ies80 = [6.536332e-9,-1.120083e-6,1.001685e-4,-9.09529e-3,6.793952e-2,999.842594];
g = 9.81;
V_m = 22.414; % (L O2 / mol O2) O2 molar volume at STP (273.15 K, 1 atm)
O2_ppm = 210000; %O2 concentration in air (mumol O2 / mol air)

% pittuus = 20;
% dz=1;
% Tz = (pittuus:-1:1)';
% zz = 0:pittuus-1;
% Pa = 1013.25;
% O2z = 9000*ones(pittuus,1);

hydr_press = NaN*ones(length(Tz),1);

density = polyval(ies80,max(0,Tz))+min(Tz,0); %kg/m^3
% Note: in equations of density it is assumed that every supercooled degree
%lowers density by 1 kg m-3 due to frazil ice formation
%(probably no practical meaning, but included for "safety")

Pa = 0.98692e-3*Pa; %atm

Tz = max(0,Tz);

% Hydrostatic pressure (kPa),calculated in the middle of the layer

hydr_press(1) = g*density(1)*dz/2;
dummy = hydr_press(1);
for i=2:length(hydr_press)
    hydr_press(i) = dummy + g*(density(i-1)+density(i))/2*dz; 
    dummy = hydr_press(i);
end

% Total pressure in the middle of the layers
tot_press = Pa + 0.98692e-5.*hydr_press; %atm

%O2 solubility constant

lnbeeta = -58.3877+85.8079*(100./(Tz+273.15))+23.8439*log((Tz+273.15)./100);
beeta = exp(lnbeeta); %Bunsen solubility coefficient (L O2 / (L water * atm))
K0_O2 = beeta./V_m; %(mol O2 / (L water * atm))

%O2 equilibrium concentration (g/mol * mol/(L*atm) * mumol/mol * atm) = mug/l = mg/m^3
  % on the surface
O2_eq_surf = 32*K0_O2*O2_ppm*Pa;
  % in the middle of the layers
O2_eq_abs = 32*K0_O2*O2_ppm.*tot_press;

% Relative oxygen concentration

% Including hydrostatic pressure
O2_sat_abs = O2z./O2_eq_abs;

% Relative to surface
O2_sat_rel = O2z./O2_eq_surf;

end