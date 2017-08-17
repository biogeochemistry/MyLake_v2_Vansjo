%% new_pH_module: estimates pH based on electro-neutrality
function [pH] = new_pH_module(acids, ions)
    pH = 7;


%% Acid_1: returns charge of sum of acid at a given pH
function [charge] = Acid_1(pH, Ct, pK)
    charge = 0

%% Acid_2: returns charge sum of diprotonic acid at a given pH
function [charge] = Acid_2(pH, Ct, pK1, pK2)
    charge = ;

