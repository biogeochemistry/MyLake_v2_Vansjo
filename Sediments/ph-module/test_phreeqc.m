clear all
mex -v CC=gcc LD=g++ LINKOPTIM='-O3' COPTIMFLAGS='-O3 -DNDEBUG -ffast-math' -I/usr/local/include/ -I/Users/MarkelovIgor/eigen -L/usr/local/lib/ -liphreeqc pH_phreeqc.cpp
% in=[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ]
% [out]=pH(in);

% clear all; close all; clc;
% tic;
format long;
n =64;
x=linspace(0,2,n);

M=4;
u0=M*-sin(pi*x);
u1=M*sin(pi*x);
m = 2;

Kc1=5.01*10^(-7).*10^6; Kc2=4.78*10^(-11).*10^6; Knh=5.62*10^(-10).*10^6; Khs=1.3*10^(-7).*10^6; Kw=10^(-14).*10^12; Kc0 = 1.7*10^(-3); 

% [H2CO3]/[CO2] ≈ 1.7×10−3

% H = load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/H.mat'); % H
H=zeros(n,m);
OH = zeros(n,m);
HS = load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/HS.mat'); % HS
HS = HS.S1;
H2S = load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/H2S.mat'); % H2S
H2S = H2S.S2;
NH3 = zeros(n,m);
load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/NH4.mat')
CO3 = zeros(n,m);
HCO3 = zeros(n,m);
CO2 = zeros(n,m);
H2CO3 = zeros(n,m);
FeS = zeros(n,m);
FeS = zeros(n,m);
FeOH3 = zeros(n,m);
FeOOH = zeros(n,m);
Ca3PO4 = zeros(n,m);
PO4adsa = zeros(n,m);
PO4adsb = zeros(n,m);

FeS(:,1) = ones(64,1)*0.001;
FeS2(:,1) = ones(64,1)*0.0015;
FeOH3(:,1) = ones(64,1)*0.0016;
FeOOH(:,1) = ones(64,1)*0.0017;
Ca3PO4(:,1)= ones(64,1)*0.0018;
PO4adsa = ones(64,1)*0.000018;
PO4adsb = ones(64,1)*0.000019;



H(:,1) = ones(64,1)*0.38;
OH(:,1) = Kw./H(:,1);
H2S(:,1) = H2S(:,1);
HS(:,1) = H2S(:,1).*Khs ./ H(:,1);
NH3(:,1) = Knh.*NH4./H(:,1);
NH4(:,1) = NH4(:,1);
CO2(:,1) = 10.^u0';
H2CO3(:,1) = 1./CO2(:,1);
HCO3(:,1) = Kc1 .* H2CO3(:,1) ./ H(:,1);
CO3(:,1) = Kc2 .* HCO3(:,1) ./ H(:,1);


load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/Fe2.mat'); % Fe2

load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/Ca2.mat'); % Ca2

load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/NO3.mat'); % NO3

load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/SO4.mat'); % SO4

load('/Users/MarkelovIgor/Documents/git/biogeochemistry/equilibria/init_concentrations/PO4.mat'); % PO4

% break
% tic
for i=1:1
    in =[H(:,1) HCO3(:,1) CO2(:,1) CO3(:,1) NH3(:,1) NH4(:,1) HS(:,1) H2S(:,1) OH(:,1) H2CO3(:,1) Fe2(:,1) Ca2(:,1) NO3(:,1) SO4(:,1) PO4(:,1) FeS(:,1) FeS2(:,1) FeOH3(:,1) FeOOH(:,1) Ca3PO4(:,1) PO4adsa(:,1) PO4adsb(:,1)];
    tic;
    [out] = pH_phreeqc(64,in)
    toc
    % H(i,2) = out(1)*out(1);
    % HCO3(i,2) = out(2)*out(2);
    % CO2(i,2) = out(3)*out(3);
    % CO3(i,2) = out(4)*out(4);
    % NH3(i,2) = out(5)*out(5);
    % NH4(i,2) = out(6)*out(6);
    % HS(i,2) = out(7)*out(7);
    % H2S(i,2) = out(8)*out(8);
    % OH(i,2) = out(9)*out(9);
    % H2CO3(i,2) = out(10)*out(10);
    % Fe2(i,2) = out(10)*out(11);
    % Ca2(i,2) = out(11)*out(12);
    % NO3(i,2) = out(12)*out(13);
    % SO4(i,2) = out(13)*out(14);
    % PO4(i,2) = out(14)*out(15);

end

% toc

% Constants
% Kc1_eq=10^(-6.4); Kc2_eq=10^(-10.3)*10^6; Knh_eq=10^(-9.3)*10^6; Khs_eq=10^(-7)*10^6; Kw_eq=10^(-14)*10^12;

% K_w = H.*OH;
% K_c0 = H2CO3 ./ CO2;
% K_c1 = H.*HCO3 ./ H2CO3;
% K_c2 = H.*CO3./HCO3;
% K_nh = H.*NH3 ./ NH4;
% K_hs = H.*HS ./ H2S;

% % charge balance
% CB = H(:,2) + NH4(:,2) + 2*Fe2 + 2*Ca2 - (HCO3(:,2) + 2*(CO3(:,2)) + HS(:,2) + OH(:,2)  + NO3 + 2*SO4 + 3*PO4 );
% CB = H + NH4  - (HCO3 + 2*(CO3) + HS + OH );


% Mol_fr_CO2 = CO2 ./ (CO3+CO2+HCO3+H2CO3);
% Mol_fr_CO3 = CO3 ./ (CO3+CO2+HCO3+H2CO3);
% Mol_fr_HCO3 = HCO3 ./ (CO3+CO2+HCO3+H2CO3);
% Mol_fr_H2CO3 = H2CO3 ./ (CO3+CO2+HCO3+H2CO3);

% Mol_fr_NH3 = NH3 ./ (NH3+NH4);
% Mol_fr_NH4 = NH4 ./ (NH3+NH4);

% Mol_fr_HS = HS ./ (HS + H2S);
% Mol_fr_H2S = H2S ./ (HS + H2S);

% pH = -log10(H(:,end)*10^-6);
% pH0 = -log10(H(:,1)*10^-6);
% figure
% semilogy(pH0,H(:,1),'x',pH0,OH(:,1),'x',pH0,H2S(:,1),'.',pH0,HS(:,1),'.', pH0,NH3(:,1),'+',pH0,NH4(:,1),'+', pH0,HCO3(:,1),'s',pH0,CO3(:,1),'s',pH0,CO2(:,1),'s')
% legend('H', 'OH', 'H2S', 'HS', 'NH3', 'NH4', 'HCO3', 'CO3', 'CO2')
% xlabel('pH')
% ylabel('Concentration in umol/L')
% title('Initial concentrations of species')
% xlim([0 14])


% figure
% semilogy(pH,H(:,end),'x',pH,OH(:,end),'x',pH,H2S(:,end),'.',pH,HS(:,end),'.', pH,NH3(:,end),'+',pH,NH4(:,end),'+', pH,HCO3(:,end),'s', pH,H2CO3(:,end),'s',pH,CO3(:,end),'s',pH,CO2(:,end),'s')
% legend('H', 'OH', 'H2S', 'HS', 'NH3', 'NH4', 'HCO3', 'H2CO3', 'CO3', 'CO2')
% xlabel('pH')
% ylabel('Concentration in umol/L')
% title('Eq concentrations of species')
% xlim([0 14])


% figure
% plot(pH,Mol_fr_CO2(:,end),'x', pH,Mol_fr_CO3(:,end),'x', pH,Mol_fr_HCO3(:,end),'x',pH,Mol_fr_H2CO3(:,end),'x', pH, Mol_fr_NH3(:,end),'o', pH,Mol_fr_NH4(:,end),'o' , pH, Mol_fr_HS(:,end),'d', pH,Mol_fr_H2S(:,end),'d' )
% ylim([0 1])
% legend('CO2','CO3','HCO3','H2CO3','NH3', 'NH4', 'HS', 'H2S')
% xlabel('pH')
% ylabel('Mole fraction')
% title('Eq Mole fractions')
% xlim([0 14])


