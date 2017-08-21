function [alphas] = alpha(pH, pKs)

    Ks = 10.^(-pKs);
    h3o = 10.^(-pH);

    if pKs == 0
        alphas = 1;
    elseif size(pKs, 2) == 1
        a1 = h3o ./ (h3o + Ks(1)*h3o);
        a2 = 1 - a1;
        alphas = [a1, a2];
    elseif size(pKs, 2) == 2
        a1 = h3o.^2./ (h3o.^2 + Ks(1)*h3o + Ks(1)*Ks(2));
        a2 = h3o.*Ks(1) ./ (h3o.^2 + Ks(1)*h3o + Ks(1)*Ks(2));
        a3 = 1 - a1 - a2;
        alphas = [a1, a2, a3];

    elseif size(pKs, 2) == 3
        a1 = h3o.^3./ (h3o.^3 + Ks(1)*h3o.^2 + Ks(1)*Ks(2)*h3o + Ks(1)*Ks(2)*Ks(3));
        a2 = Ks(1)*h3o.^2./ (h3o.^3 + Ks(1)*h3o.^2 + Ks(1)*Ks(2)*h3o + Ks(1)*Ks(2)*Ks(3));
        a3 = Ks(2)*Ks(1)*h3o./ (h3o.^3 + Ks(1)*h3o.^2 + Ks(1)*Ks(2)*h3o + Ks(1)*Ks(2)*Ks(3));
        a4 = 1 - a1 - a2 - a3;
        alphas = [a1, a2, a3, a4];
    else
        % # These are the powers that the H3O+ concentrations will be raised.
        powers = [0:size(pKs, 2)];
        powers_rev = powers(end:-1:1);

         % # Calculate the H3O+ concentrations raised to the powers calculated
        % # above (in reverse order).
        h3o_pow = bsxfun(@power, h3o, powers_rev);

        % # Calculate a cumulative product of the Ka values. The first value
        % # must be 1.0
        Ka_prod = cumprod([1,Ks]);
        % # Multiply the H3O**power values times the cumulative Ka product.
        h3o_Ka = bsxfun(@times,h3o_pow, Ka_prod);
        den = sum(h3o_Ka, 2);
        alphas =bsxfun(@rdivide, h3o_Ka, den);
    end
end
