function [CO2aq, HCO3aq, CO3aq, CO2frac, HCO3frac, CO3frac] = carbonequilibrium(DICz,Tz,pH)
%Calculates the fractions of dissolved inorganic carbon
%Inputs
%pH          pH
%DICzt       DICzt profile (mg/m^3)
%Tz          water temperature profile (deg C)



%Outputs
%CO2aq       carbon dioxide (mg/m^3)
%CO2frac     CO2 fraction (mol CO2 / mol DIC)
%pCO2        CO2 partial pressure (atm)
%HCO3frac    HCO3 fraction (mol HCO3 / mol DIC)
%CO3frac     CO3 fraction (mol CO3 / mol DIC)
%alkalinity  alkalinity (mg/l)
%C_acid      acid carbon (mg/m^3)
%C_basic     basic carbon (mg/m^3)

%Water density polynomial

global ies80;

density = (polyval(ies80,max(0,Tz))+min(Tz,0))*0.001; %kg/l
% Note: in equations of density it is assumed that every supercooled degree lowers density by
% 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")

DICz = 0.001.*DICz; %mg/l
Tz = max(0,Tz)+273.15;

%Carbon acid solubility and dissociation constants (Millero, 1995)

K0 = -60.2409+93.4517.*(100./Tz)+23.3585*log(Tz/100); %~mol/(kg*atm)
K0 = exp(K0).*density; %CO2; mol/(kg*atm)*kg/l = mol/(l*atm)
K1 = 290.9097-14554.21./Tz-45.0575*log(Tz); %~mol/kg
K1 = exp(K1); %HCO3; mol/kg
K2 = 207.6548-11843.79./Tz-33.6485*log(Tz); %~mol/kg
K2 = exp(K2); %CO3; mol/kg
Kw = 148.9802-13847.26./Tz-23.6521*log(Tz); %~mol/kg
Kw = exp(Kw); %H2O; (mol/kg)^2

%Hydrogen ion concentration

Hplus = NaN*zeros(length(pH),1);
for i=1:length(pH)
Hplus(i) = 10^(-pH(i)); %mol/l
end
Hplus = Hplus./density; %mol/kg

%fraction of carbonates in a particular form
%http://www.chem.usu.edu/~sbialkow/Classes/3650/Carbonate/Carbonic%20Acid.html

CO2mfrac = Hplus.*Hplus./((Hplus.*Hplus+Hplus.*K1+K1.*K2)); %mol CO2 / mol DIC
HCO3mfrac = Hplus.*K1./((Hplus.*Hplus+Hplus.*K1+K1.*K2)); %-mol HCO3 / mol DIC
CO3mfrac = K1.*K2./((Hplus.*Hplus+Hplus.*K1+K1.*K2)); %-mol CO3 / mol DIC

M_DIC = CO2mfrac*44.01+HCO3mfrac*61.01+CO3mfrac*60.01; % DIC molar mass g/mol

DICmz = DICz./M_DIC; %DIC concentration in mmol DIC/l

CO2frac = CO2mfrac*44.01./M_DIC; % mg CO2 / mg DIC
HCO3frac = HCO3mfrac*61.01./M_DIC; % mg HCO3 / mg DIC
CO3frac = CO3mfrac*60.01./M_DIC; % mg CO3 / mg DIC

%Carbon dioxide partial pressure & concentration mg/m^3

pCO2 = CO2mfrac.*DICmz./(1000.*K0); %mol CO2/mol DIC * mmol DIC/l * l*atm/mmol CO2 = atm
CO2aq = 1000*1000*44.01.*K0.*pCO2; %mug/mol * mol/l = mug/l = mg/m^3
HCO3aq = 1000*HCO3frac.*DICz; %mug/l = mg/m^3
CO3aq = 1000*CO3frac.*DICz; %mug/l = mg/m^3

%Alkalinity mmol/l

%OH = Kw./Hplus.*density; %mol/l
%Hpluss = Hplus.*density; %mol/l

%alkalinitym = (HCO3mfrac + 2*CO3mfrac).*DICmz+OH*1000-Hpluss*1000; %mmol/l

%Alkalinity mg/l
%alkalinity = (HCO3frac + 2*CO3frac).*DICz+OH*1000*17-Hpluss*1000*1.008; %mg/l

%C_acid = CO2aq+0.5*HCO3aq+0.5*(1e6*Hpluss*1.008-1e6*OH*17); %mg/m^3
%C_basic = CO3aq+0.5*HCO3aq-0.5*(1e6*Hpluss*1.008-1e6*OH*17); %mg/m^3

end

