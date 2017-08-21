%% new_pH_module: estimates pH based on electro-neutrality
function [x] = new_pH_module(Aquatic_system, pHs)

    x = estime_charge(pHs);
    [M,I] = min(x);
    x = pHs(I);

    function [ch] = estime_charge(pH)
        x = 10.^(-pH) -  (10.^(-14+pH));
        fields = fieldnames(Aquatic_system);
        for i = 1:numel(fields)
            conc = Aquatic_system.(fields{i}).conc;
            charge = Aquatic_system.(fields{i}).charge;
            alphas = alpha(pH, Aquatic_system.(fields{i}).pKs);
            x = x + sum(bsxfun(@times,bsxfun(@times,conc, charge), alphas), 2);
            % x = x + sum(conc .* alphas .* charge);
          % x += sum(Aquatic_system.(fields{i}).conc .* Aquatic_system.(fields{i}).charge .* alpha(pH, Aquatic_system.(fields{i}).pKs)
        end
        ch = abs(x);
    end

end

function [x] = test()
    acid1 = acid([3.6, 10.32], 0, 0.0000001);

    aq_sys.ca = neutral(2, 1e-6);

    aq_sys.acid1 = acid1;
    aq_sys.acid2 = acid1;


    pHs = linspace(0,14,1400)';
    x = new_pH_module(aq_sys, pHs);
end






