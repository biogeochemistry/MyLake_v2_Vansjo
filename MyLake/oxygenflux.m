function [O2_new,surfaceflux,O2_eq,K0_O2] = oxygenflux(O2_old,Ws,Pa,T0,dz)
%Oxygen surface flux

%Inputs
%O2_old      Old O2 concentration in the surface layer (mg/m^3) 
%Ws          Wind speed (m/s)
%Pa          Air pressure (mbar = hPa)  
%T0          Water surface temperature
%dz          Surface layer depth

%Outputs
%O2_new      New O2 concentration (mg/m^3)
%surfaceflux O2 air-water surface flow mg/(m^2*d)
%O2_eq       Oxygen equilibrium concentration (mg/m^3)
%K0_O2       Oxygen solubility constant (mol O2 / (L water * atm))

V_m = 22.414; % (L O2 / mol O2) O2 molar volume at STP (273.15 K, 1 atm)

Pa = 0.98692e-3*Pa; %atm
O2_ppm = 210000; %O2 concentration in air (mumol O2 / mol air)

if(T0 < 0)
    T0 = 0;
end    

%O2 solubility constant

lnbeeta = -58.3877+85.8079*(100./(T0+273.15))+23.8439*log((T0+273.15)./100);
beeta = exp(lnbeeta); %Bunsen solubility coefficient (L O2 / (L water * atm))
K0_O2 = beeta./V_m; %(mol O2 / (L water * atm))

%O2 equilibrium concentration (g/mol * mol/(L*atm) * mumol/mol * atm) = mug/l = mg/m^3 

O2_eq = 32*K0_O2*O2_ppm*Pa;

alpha = 1; %Chemical enhancement factor
A = 1800.6; %Schmidt number polynomial fit coefficients
B = 120.10;
C = 3.7818;
D = 0.047608;

k_600 = 2.07+0.215*Ws^1.7; %transfer velocity for Schmidt number 600 (cm/h)
schmidt = A-B*T0+C*T0.*T0-D*T0.*T0.*T0; %Schmidt number for O2 (-)
k_O2 = k_600*(schmidt/600)^(-0.666); %transfer velocity for CO_2 (cm/h)

surfaceflux = 0.24*alpha*k_O2.*(O2_old-O2_eq); %(m/cm * h/d * cm/h *mg/m^3 = mg/(m^2*d)) 
O2_new = O2_old-surfaceflux/dz; %(mg/m^3); time step = 1 d)

%alpha = beeta.*(1 + 0.00367.*T0);
% K0_Trumans = (0.046*(T0+273.15).*(T0+273.15) + 203.357*(T0+273.15).*log((T0+273.15)./298) -...
%              (299.378 + 0.092*(T0+273.15)).*((T0+273.15)-298) - 20.591e3) ./ (8.3144*(T0+273.15));
% K0_Trumans =exp(K0_Trumans);         
         
%Bunsenin saa nyt K_0:ksi ja sen edelleen tiheyden kanssa
%skaalattua oikeaan yksikköön. Suora Trumansin kaava antaa saman tuloksen
%(molality eli mol/kg!) kuin Bunsenin kautta puljattu.
         
%http://www.medilexicon.com/medicaldictionary.php?t=18715
%Bunsenin ja Ostwaldin oikea (?) suhde.