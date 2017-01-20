function [CO2_new,surfaceflux,CO2_eq,K0,CO2_ppm] = carbondioxideflux(CO2_old,Ws,Pa,T0,dz,date)
%Carbon dioxide surface flux

%Inputs
%CO2_old     Old CO2 concentration in the surface layer (mg/m^3)
%Ws          Wind speed (m/s)
%Pa          Air pressure (mbar = hPa)
%T0          Water surface temperature
%dz          Surface layer depth
%date        Date number

%Outputs
%CO2_new     New CO2 concentration (mg/m^3)
%surfaceflux CO2 air-water surface flux mg/(m^2*d)
%K0          Carbon dioxide solubility coefficient (mol/(l*atm))
%CO2_ppm     Carbon dioxide concentration in air (mumol CO2 / mol air)

%Water density polynomial

global ies80;

density = (polyval(ies80,max(0,T0))+min(T0,0))*0.001; %kg/l
% Note: in equations of density it is assumed that every supercooled degree lowers density by
% 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")

Pa = 0.98692e-3*Pa; %atm
T0 = max(T0,0);

%CO2 solubility constant

K0 = -60.2409+93.4517.*(100./(T0+273.15))+23.3585*log((T0+273.15)/100); %~mol/(kg*atm)
K0 = exp(K0)*density; %CO_2; mol/(kg*atm)*kg/l = mol/(l*atm)

%CO2 concentration in air (mumol CO2 / mol air)
% TODO: This returns vector and not a single value.
CO2_ppm = (380+4*sin(((date+0.310625)/365.2425*2*pi)));

%CO2 equilibrium concentration (g/mol*mol/(l*atm)*(mumol/mol)*atm = mug/l = mg/m^3)

CO2_eq = 44.01*K0*CO2_ppm*Pa;

alpha = 1; %Chemical enhancement factor
A = 1911.1; %Schmidt number polynomial fit coefficients
B = 118.11;
C = 3.4527;
D = 0.041320;

k_600 = 2.07+0.215*Ws^1.7; %transfer velocity for Schmidt number 600 (cm/h)
schmidt = A-B*T0+C*T0^2-D*T0^3; %Schmidt number for CO2 (-)
k_CO2 = k_600*(schmidt/600)^(-0.5); %transfer velocity for CO2 (cm/h)

surfaceflux = 0.24*alpha*k_CO2*(CO2_old-CO2_eq); %(m/cm * h/d * cm/h *mg/m^3 = mg/(m^2*d))
CO2_new = CO2_old-surfaceflux/dz; %(mg/m^3); time step = 1 d)

end
