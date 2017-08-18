function [acid] = acid(pKs, charge, conc)
    acid.pKs = pKs;
    acid.charge = [charge:-1:charge - size(pKs,2)];
    acid.conc = conc;
end
