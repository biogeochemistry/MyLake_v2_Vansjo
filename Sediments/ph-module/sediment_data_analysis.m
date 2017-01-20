% clear all;
% load('sediments_data_basin3.mat')

HCO3=sediments_data_basin3{28, 1};
CO2=sediments_data_basin3{26, 1};
CO3=sediments_data_basin3{27, 1};
NH3=sediments_data_basin3{29, 1};
NH4=sediments_data_basin3{20, 1};
H2CO3= sediments_data_basin3{30, 1};

H2S =(sediments_data_basin3{6, 1});
HS =(sediments_data_basin3{7, 1});
H =(sediments_data_basin3{21, 1});
OH=(sediments_data_basin3{25, 1});
Fe2=(sediments_data_basin3{5, 1});
Ca2=(sediments_data_basin3{22, 1});
NO3=(sediments_data_basin3{19, 1});
SO4=(sediments_data_basin3{4, 1});
PO4=(sediments_data_basin3{16, 1});

K_w = H.*OH;
K_c0 = H2CO3 ./ CO2;
K_c1 = H.*HCO3 ./ H2CO3;
K_c2 = H.*CO3./HCO3;
K_nh = H.*NH3 ./ NH4;
K_hs = H.*HS ./ H2S;
Kc1=5.01*10^(-7).*10^6; Kc2=4.78*10^(-11).*10^6; Knh=5.62*10^(-10).*10^6; Khs=1.3*10^(-7).*10^6; Kw=10^(-14).*10^12; Kc0 = 1.7*10^(-3); 

% charge balance
CB = H + NH4 + 2*Fe2 + 2*Ca2 - (HCO3 + 2*(CO3) + HS + OH  + NO3 + 2*SO4 + 3*PO4 );
% CB = H + NH4  - (HCO3 + 2*(CO3) + HS + OH );


Mol_fr_CO2 = CO2 ./ (CO3+CO2+HCO3+H2CO3);
Mol_fr_CO3 = CO3 ./ (CO3+CO2+HCO3+H2CO3);
Mol_fr_HCO3 = HCO3 ./ (CO3+CO2+HCO3+H2CO3);
Mol_fr_H2CO3 = H2CO3 ./ (CO3+CO2+HCO3+H2CO3);

Mol_fr_NH3 = NH3 ./ (NH3+NH4);
Mol_fr_NH4 = NH4 ./ (NH3+NH4);

Mol_fr_HS = HS ./ (HS + H2S);
Mol_fr_H2S = H2S ./ (HS + H2S);

pH = -log10(H*10^-6);

% figure
% semilogy(pH,H,'x',pH,OH,'x',pH,H2S,'.',pH,HS,'.', pH,NH3,'+',pH,NH4,'+', pH,HCO3,'s',pH,CO3,'s',pH,CO2,'s')
% legend('H', 'OH', 'H2S', 'HS', 'NH3', 'NH4', 'HCO3', 'CO3', 'CO2')
% xlabel('pH')
% ylabel('Concentration in umol/L')
% title('Equilibrium concentrations of species')
% xlim([0 14])


% figure
% plot(pH,Mol_fr_CO2,'x', pH,Mol_fr_CO3,'x', pH,Mol_fr_HCO3,'x', pH, Mol_fr_NH3,'o', pH,Mol_fr_NH4,'o' , pH, Mol_fr_HS,'d', pH,Mol_fr_H2S,'d' )
% ylim([0 1])
% legend('CO2','CO3','HCO3','NH3', 'NH4', 'HS', 'H2S')
% xlabel('pH')
% ylabel('Mole fraction')
% title('Equilibrium Mole fractions')
% xlim([0 14])