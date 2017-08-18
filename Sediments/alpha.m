function [alphas] = alpha(pH, pKs)

    if pKs == 0
        alphas = 1;
    else
        h3o = 10.^(-pH);
        % # These are the powers that the H3O+ concentrations will be raised.
        powers = [0:size(pKs, 2)];
        powers_rev = powers(end:-1:1);

         % # Calculate the H3O+ concentrations raised to the powers calculated
        % # above (in reverse order).
        h3o_pow = bsxfun(@power, h3o, powers_rev);
        Kas = 10.^(-pKs);
        % # Calculate a cumulative product of the Ka values. The first value
        % # must be 1.0
        Ka_prod = cumprod([1,Kas]);
        % # Multiply the H3O**power values times the cumulative Ka product.
        h3o_Ka = bsxfun(@times,h3o_pow, Ka_prod);
        den = sum(h3o_Ka, 2);
        alphas =bsxfun(@rdivide, h3o_Ka, den);
    end
end
